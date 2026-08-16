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

  ## Psi_s is always diagonal
  inv_Psi_s <- lapply( Psi_s, .inv_Psi_vec )

  ## delta_s = (I_k + Phi' Psi_s^-1 Phi)^-1 Phi' Psi_s^-1
  ## a single Woodbury application per study
  delta_Phi <- Map( .wb_identity, list( Phi ), inv_Psi_s, list( I_k ) )

  Txsfcs     <- Map( .get_exp_xf, delta_Phi, cov_s )
  Tfcsfcs    <- Map( .get_exp_ff, list( Phi ), delta_Phi, cov_s, list( I_k ) )
  E_fis_x_is <- Map( .get_exp_f, delta_Phi, X_s_tilde )

  out <- list( Txsfcs = Txsfcs, Tfcsfcs = Tfcsfcs, E_fis_x_is = E_fis_x_is )

  if ( getdet ) {
    ## Sig_s^-1 = Psi_s^-1 - Psi_s^-1 Phi delta_s   (dense p x p, via Woodbury)
    out$Sig_s1 <- Map( .wb_identity2, inv_Psi_s, list( Phi ), list( I_k ) )

    ## det(Sig_s) via the matrix-determinant lemma: det(Psi_s) * det(I_k + Phi' Psi_s^-1 Phi)
    det_Psi_s <- lapply( Psi_s, function( Psi ) prod( diag( Psi ) ) )
    inner     <- lapply( inv_Psi_s, function( iP ) I_k + crossprod( Phi * iP, Phi ) )
    out$ds_s  <- Map( function( dPsi, inn ) dPsi * det( inn ), det_Psi_s, inner )
  }

  out
}
