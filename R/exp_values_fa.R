## =============================================================================
## .exp_values_fa()  --  E-step for the study-specific FA module (ecm_fa())
##
## Plain per-study factor analysis has no common factor: x_is = Lambda_s l_is +
## e_is, e_is ~ N(0, Psi_s), l_is ~ N(0, I). Unlike ecm_msfr()/ecm_fr(), there
## is no Phi term at all, so this does *not* reuse the generic .exp_values()
## in exp_values.R (which always carries a Phi block and a double Woodbury
## reduction for it) -- nor the "Phi = 0 dummy column" trick ecm_fa() used to
## rely on to borrow that generic function. That trick still pays for a
## Woodbury reduction against a rank-0 factor on every call, and it broke
## outright once .exp_values() was refactored for MSFR (ecm_fa.R's call,
## `.exp_values(Phi, Omega_s, Psi_s, Psi_s1, cov_s, getdet = TRUE)`, no longer
## matches that function's signature or return value).
##
## Only a single Woodbury application per study is needed here (versus the
## double reduction .exp_values() performs to account for Phi), and only
## Txsfs/Tfsfs are returned (plus Sig_s1/ds_s when getdet = TRUE): these are
## the only fields ecm_fa() actually reads from its E-step call. The generic
## .exp_values() also returns Txsfcs/Tfcsfcs/Tfcsfs -- cross terms involving
## Phi -- which ecm_fa() unpacks but never uses (dead reads, a leftover of the
## Phi = 0 trick); they have no meaning in a Phi-free model and are dropped.
##
## Sig_s1/ds_s are returned as a dense p x p matrix / a plain determinant
## (not the compact Woodbury representation used in the newer
## exp_values_fr.R), because ecm_fa()'s stopping rule calls the *existing*
## .loglik_ecm(Sig_s1, ds_s, n_s, cov_s), which does Sig_s1[[s]] %*% cov_s[[s]]
## and log(ds_s[[s]]) directly -- it has no notion of a compact
## representation, and changing that is out of scope here. det(Sig_s) still
## avoids an O(p^3) dense determinant via the matrix-determinant lemma
## (Sig_s is p x p, but Psi_s is diagonal and Lambda_s is p x q with q << p):
##   det(Psi_s + Lambda_s Lambda_s') = det(Psi_s) * det(I_q + Lambda_s' Psi_s^-1 Lambda_s)
## =============================================================================

#' Expected values for the FA E-step (used by \code{ecm_fa()})
#'
#' Computes the sufficient statistics needed by the CM-steps of \code{ecm_fa()}
#' for a set of independent, per-study factor analysis models
#' \eqn{x_{is} = \Lambda_s l_{is} + e_{is}}{x_is = Lambda_s l_is + e_is}, with
#' \eqn{e_{is} \sim N(0, \Psi_s)}{e_is ~ N(0, Psi_s)} and \eqn{l_{is} \sim N(0, I)}{l_is ~ N(0, I)}.
#' There is no common factor in this model, so unlike \code{.exp_values()} (used by
#' \code{ecm_msfr()}), no \code{Phi} argument is needed or computed.
#' @param Lambda_s List of length \eqn{S}{S}, study-specific loading matrices, each
#' \eqn{p \times q_s}{p x q_s}.
#' @param Psi_s List of length \eqn{S}{S}, study-specific idiosyncratic covariance matrices, each a
#' \eqn{p \times p}{p x p} diagonal matrix.
#' @param cov_s List of length \eqn{S}{S}, study-specific sample covariance (or correlation) matrices.
#' @param getdet If \code{TRUE}, also returns \code{Sig_s1} (\eqn{\Sigma_s^{-1}}{Sigma_s^-1}) and
#' \code{ds_s} (\eqn{\det(\Sigma_s)}{det(Sigma_s)}), as required by \code{.loglik_ecm()}. Default is
#' \code{FALSE}: these are only needed for the log-likelihood/stopping-rule evaluation, not for the
#' CM parameter updates, so skipping them saves a Woodbury pass on every CM-step call.
#' @return A list containing:
#' \item{Txsfs}{list of \eqn{E[x_s l_s']}{E[x_s l_s']}, i.e. \eqn{p \times q_s}{p x q_s} matrices.}
#' \item{Tfsfs}{list of \eqn{E[l_s l_s']}{E[l_s l_s']}, i.e. \eqn{q_s \times q_s}{q_s x q_s} matrices.}
#' \item{Sig_s1}{(only if \code{getdet = TRUE}) list of \eqn{\Sigma_s^{-1}}{Sigma_s^-1}, dense
#' \eqn{p \times p}{p x p} matrices.}
#' \item{ds_s}{(only if \code{getdet = TRUE}) list of \eqn{\det(\Sigma_s)}{det(Sigma_s)} scalars.}
#' @details
#' Reuses \code{.inv_Psi_vec()}, \code{.wb_identity()} and \code{.wb_identity2()} (\code{helpers.R} /
#' \code{exp_values.R}) to apply the Woodbury identity, and \code{.get_exp_xl()} /
#' \code{.get_exp_ll()} (\code{exp_values.R}) to assemble the sufficient statistics, so none of
#' this linear algebra is duplicated across the FA/FR/MSFR modules. \code{Psi_s^{-1}} is passed to
#' \code{.wb_identity()}/\code{.wb_identity2()} as a plain vector rather than a dense matrix, since
#' it's always diagonal -- this lets them take their row-scaling fast path instead of a dense
#' \eqn{p \times p}{p x p} multiply.
#' @references De Vito, R., Bellio, R., Trippa, L. and Parmigiani, G. (2019). Multi-study Factor
#' Analysis. Biometrics, 75, 337-346.
#' @keywords internal
.exp_values_fa <- function( Lambda_s, Psi_s, cov_s, getdet = FALSE )
{
  I_j <- lapply( Lambda_s, function( L ) diag( 1, ncol( L ) ) )

  ## Psi_s is always diagonal, so its inverse is kept as a length-p vector
  ## rather than a dense p x p matrix: .wb_identity()/.wb_identity2() take
  ## the vector fast path for it (see their docs), turning what would be an
  ## O(p^2 q) dense multiply into an O(p q) row-scaling one.
  inv_Psi_s <- lapply( Psi_s, .inv_Psi_vec )

  ## delta_s = (I_qs + Lambda_s' Psi_s^-1 Lambda_s)^-1 Lambda_s' Psi_s^-1
  ## a single Woodbury application per study -- no Phi to reduce against, so
  ## unlike the generic .exp_values() there is no wb_f/wb_l double reduction.
  delta_Lambda <- Map( .wb_identity, Lambda_s, inv_Psi_s, I_j )

  Txsfs <- Map( .get_exp_xl, delta_Lambda, cov_s )
  Tfsfs <- Map( .get_exp_ll, Lambda_s, delta_Lambda, cov_s, I_j )

  out <- list( Txsfs = Txsfs, Tfsfs = Tfsfs )

  if ( getdet ) {
    ## Sig_s^-1 = Psi_s^-1 - Psi_s^-1 Lambda_s delta_s   (dense p x p, via Woodbury)
    out$Sig_s1 <- Map( .wb_identity2, inv_Psi_s, Lambda_s, I_j )

    ## det(Sig_s) via the matrix-determinant lemma: det(Psi_s) * det(I_q + Lambda_s' Psi_s^-1 Lambda_s),
    ## i.e. an O(p) diagonal product and a determinant on a q x q matrix, instead of an O(p^3)
    ## dense determinant on Sig_s itself. iP is a vector here (see above), so the inner
    ## quadratic form is assembled via row-scaling (L * iP) instead of a dense p x p multiply.
    det_Psi_s <- lapply( Psi_s, function( Psi ) prod( diag( Psi ) ) )
    inner     <- Map( function( L, iP, I ) I + crossprod( L * iP, L ), Lambda_s, inv_Psi_s, I_j )
    out$ds_s  <- Map( function( dPsi, inn ) dPsi * det( inn ), det_Psi_s, inner )
  }

  out
}
