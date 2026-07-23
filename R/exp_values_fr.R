## =============================================================================
## .exp_values_fr()  --  E-step for the factor-regression module (ecm_fr())
##
## FR has a single shared common factor and no per-study loadings:
## x_is = beta b_is + Phi f_is + e_is, e_is ~ N(0, Psi_s), f_is ~ N(0, I).
## (x_is here is already the covariate-adjusted x_tilde = x - beta*b, computed
## by ecm_fr() before calling this function.)
##
## This used to carry its own bespoke Woodbury implementation
## (.get_diag_vec/.woodbury_block/.apply_Sig_s1) and its own compact
## representation of Sig_s1, all local to this file. None of that is needed:
## the same Woodbury reduction is already implemented once, generically, in
## exp_values.R/helpers.R (.inv_Psi, .wb_identity, .wb_identity2) and reused
## as-is by .exp_values_fa() for the FA module. FR needs exactly one Woodbury
## application per study to reduce against Phi (like FA reduces against
## Lambda_s) -- there's no second factor block to double-reduce against, so
## none of .exp_values()'s Phi/Lambda_s double-reduction machinery applies
## either. .get_exp_xf()/.get_exp_ff()/.get_exp_f() (exp_values.R) already
## compute exactly the sufficient statistics FR needs from delta_Phi, so
## nothing here is re-derived that isn't already shared.
##
## Sig_s1/ds_s are dense p x p matrices / plain determinants (not the compact
## Woodbury representation the previous version of this file used), because
## ecm_fr()'s stopping rule calls the *existing* .loglik_ecm(Sig_s1, ds_s,
## n_s, cov_s), which does Sig_s1[[s]] %*% cov_s[[s]] and log(ds_s[[s]])
## directly -- same reasoning as .exp_values_fa(). det(Sig_s) still avoids an
## O(p^3) dense determinant via the matrix-determinant lemma:
##   det(Psi_s + Phi Phi') = det(Psi_s) * det(I_k + Phi' Psi_s^-1 Phi)
## =============================================================================

#' Expected values for the FR E-step (used by \code{ecm_fr()})
#'
#' Computes the sufficient statistics needed by the CM-steps of \code{ecm_fr()} for a
#' factor-regression model with a single common factor and no per-study loadings:
#' \eqn{\tilde{x}_{is} = \Phi f_{is} + e_{is}}{x_is = Phi f_is + e_is}, with
#' \eqn{e_{is} \sim N(0, \Psi_s)}{e_is ~ N(0, Psi_s)} and \eqn{f_{is} \sim N(0, I)}{f_is ~ N(0, I)},
#' where \eqn{\tilde{x}_{is}}{x_is} is the covariate-adjusted response.
#' @param Phi Common factor loading matrix, \eqn{p \times k}{p x k} (shared across all studies).
#' @param Psi_s List of length \eqn{S}{S}, study-specific idiosyncratic covariance matrices, each a
#' \eqn{p \times p}{p x p} diagonal matrix.
#' @param cov_s List of length \eqn{S}{S}, study-specific sample covariance matrices of the
#' covariate-adjusted response.
#' @param X_s_tilde List of length \eqn{S}{S}, the covariate-adjusted response data matrices
#' (\eqn{n_s \times p}{n_s x p}), used for the factor-score expectations needed by the \eqn{\beta}
#' CM-step.
#' @param getdet If \code{TRUE}, also returns \code{Sig_s1} (\eqn{\Sigma_s^{-1}}{Sigma_s^-1}) and
#' \code{ds_s} (\eqn{\det(\Sigma_s)}{det(Sigma_s)}), as required by \code{.loglik_ecm()}. Default is
#' \code{FALSE}, since these are only needed for the log-likelihood/stopping-rule evaluation.
#' @return A list containing:
#' \item{Txsfcs}{list of \eqn{E[\tilde{x}_s f_s']}{E[x_s f_s']}, i.e. \eqn{p \times k}{p x k} matrices.}
#' \item{Tfcsfcs}{list of \eqn{E[f_s f_s']}{E[f_s f_s']}, i.e. \eqn{k \times k}{k x k} matrices.}
#' \item{E_fis_x_is}{list of \eqn{E[f_{is} | x_{is}]}{E[f_is | x_is]}, i.e. \eqn{k \times n_s}{k x n_s}
#' matrices (one column per individual), used to update \eqn{\beta}.}
#' \item{Sig_s1}{(only if \code{getdet = TRUE}) list of \eqn{\Sigma_s^{-1}}{Sigma_s^-1}, dense
#' \eqn{p \times p}{p x p} matrices.}
#' \item{ds_s}{(only if \code{getdet = TRUE}) list of \eqn{\det(\Sigma_s)}{det(Sigma_s)} scalars.}
#' @details
#' Reuses \code{.inv_Psi_vec()}, \code{.wb_identity()} and \code{.wb_identity2()} (\code{helpers.R} /
#' \code{exp_values.R}) to apply the Woodbury identity, and \code{.get_exp_xf()} /
#' \code{.get_exp_ff()} / \code{.get_exp_f()} (\code{exp_values.R}) to assemble the sufficient
#' statistics -- the same building blocks \code{.exp_values_fa()} uses for the FA module.
#' \code{Psi_s^{-1}} is passed as a plain vector rather than a dense matrix, since it's always
#' diagonal -- see \code{.wb_identity()}'s docs.
#' @references Avalos-Pacheco, A., Rossell, D. and Savage, R.S. (2022). Heterogeneous Large Datasets
#' Integration Using Bayesian Factor Regression. Bayesian Analysis, 17, 33-66.
#' @keywords internal
.exp_values_fr <- function( Phi, Psi_s, cov_s, X_s_tilde, getdet = FALSE )
{
  k   <- ncol( Phi )
  I_k <- diag( 1, k )

  ## Psi_s is always diagonal, so its inverse is kept as a length-p vector
  ## rather than a dense p x p matrix: .wb_identity()/.wb_identity2() take
  ## the vector fast path for it (see their docs), turning what would be an
  ## O(p^2 k) dense multiply into an O(p k) row-scaling one.
  inv_Psi_s <- lapply( Psi_s, .inv_Psi_vec )

  ## delta_s = (I_k + Phi' Psi_s^-1 Phi)^-1 Phi' Psi_s^-1  -- a single
  ## Woodbury application per study, exactly as .exp_values_fa() does for
  ## Lambda_s; FR has no second loading matrix to reduce against.
  delta_Phi <- Map( .wb_identity, list( Phi ), inv_Psi_s, list( I_k ) )

  Txsfcs     <- Map( .get_exp_xf, delta_Phi, cov_s )
  Tfcsfcs    <- Map( .get_exp_ff, list( Phi ), delta_Phi, cov_s, list( I_k ) )
  E_fis_x_is <- Map( .get_exp_f, delta_Phi, X_s_tilde )

  out <- list( Txsfcs = Txsfcs, Tfcsfcs = Tfcsfcs, E_fis_x_is = E_fis_x_is )

  if ( getdet ) {
    ## Sig_s^-1 = Psi_s^-1 - Psi_s^-1 Phi delta_s   (dense p x p, via Woodbury)
    out$Sig_s1 <- Map( .wb_identity2, inv_Psi_s, list( Phi ), list( I_k ) )

    ## det(Sig_s) via the matrix-determinant lemma: det(Psi_s) * det(I_k + Phi' Psi_s^-1 Phi),
    ## i.e. an O(p) diagonal product and a determinant on a k x k matrix, instead of an O(p^3)
    ## dense determinant on Sig_s itself. iP is a vector here (see above), so the inner
    ## quadratic form is assembled via row-scaling (Phi * iP) instead of a dense p x p multiply.
    det_Psi_s <- lapply( Psi_s, function( Psi ) prod( diag( Psi ) ) )
    inner     <- lapply( inv_Psi_s, function( iP ) I_k + crossprod( Phi * iP, Phi ) )
    out$ds_s  <- Map( function( dPsi, inn ) dPsi * det( inn ), det_Psi_s, inner )
  }

  out
}
