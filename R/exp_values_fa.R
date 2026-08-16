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

  ## Psi_s is always diagonal
  inv_Psi_s <- lapply( Psi_s, .inv_Psi_vec )

  ## delta_s = (I_qs + Lambda_s' Psi_s^-1 Lambda_s)^-1 Lambda_s' Psi_s^-1
  ## a single Woodbury application per study
  delta_Lambda <- Map( .wb_identity, Lambda_s, inv_Psi_s, I_j )

  Txsfs <- Map( .get_exp_xl, delta_Lambda, cov_s )
  Tfsfs <- Map( .get_exp_ll, Lambda_s, delta_Lambda, cov_s, I_j )

  out <- list( Txsfs = Txsfs, Tfsfs = Tfsfs )

  if ( getdet ) {
    ## Sig_s^-1 = Psi_s^-1 - Psi_s^-1 Lambda_s delta_s   (dense p x p, via Woodbury)
    out$Sig_s1 <- Map( .wb_identity2, inv_Psi_s, Lambda_s, I_j )

    ## det(Sig_s) via the matrix-determinant lemma: det(Psi_s) * det(I_q + Lambda_s' Psi_s^-1 Lambda_s),
    det_Psi_s <- lapply( Psi_s, function( Psi ) prod( diag( Psi ) ) )
    inner     <- Map( function( L, iP, I ) I + crossprod( L * iP, L ), Lambda_s, inv_Psi_s, I_j )
    out$ds_s  <- Map( function( dPsi, inn ) dPsi * det( inn ), det_Psi_s, inner )
  }

  out
}
