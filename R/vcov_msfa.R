
#' Variance matrix of MLE estimates for a MSFA model
#'
#' Computes the inverse observed information for a MSFA model
#'
#'
#' Numerical differentiation is employed to obtain the observed information matrix at a
#' given parameter values, so that when the parameter values equals the MLE the function
#' returns the estimated variance matrix of the fitted model. The method is rather inefficient, and
#' it may lead to long computations, though the function is designed to be called only once after the
#' estimation has been carried out. However, it would be relatively straightforward to employ analytical
#' differentiation at least for the log-likelihood gradient, and this may be implemented in future
#' releases of the code.
#'
#' @param X_s List of lenght \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all
#' the studies.
#' @param mle The object returned by \code{ecm_msfa}.
#' @param getgrad Should the function return also the gradient at \code{mle}? Default is \code{FALSE}.
#' @return A list with exactly the same structure of the three slots \code{Phi}, \code{Lambda_s} and
#' \code{Psi_s} of \code{mle}, but containing the standard errors rather than the point estimates.
#' Furthemore, slots for the hessian matrix and the gradient at \code{mle} are included, the latter
#' is not NULL when \code{getgrad}  is \code{TRUE}.
#' @export
#' @import statmod
#' @importFrom pracma grad
#' @importFrom pracma hessian
vcov_msfa <- function(X_s, mle, getgrad = TRUE)
{
  constraint <- mle$constraint
  p <- ncol(X_s[[1]])
  k <- ncol(mle$Phi)
  S <- length(X_s)
  j_s <- c()
  for(s in 1:S) j_s[[s]] <- ncol(mle$Lambda_s[[s]])
  theta <- .param2vect(mle, mle$constraint)
  gout <- NULL
  if(getgrad) gout <- pracma::grad(.loglik_int, x = theta, n_s = mle$n_s,
                                  cov_s = mle$cov_s, k = k, j_s = j_s, constraint=constraint)
  hout <- pracma::hessian(.loglik_int, x = theta, n_s = mle$n_s, cov_s = mle$cov_s, k = k, j_s = j_s,
                          constraint=constraint)
  seout <- sqrt(diag(chol2inv(chol(-hout))))
  ###### re-arrange the ses like in param
  param <- .vect2param(seout, mle, constraint, p, k, j_s)
  out <- list(Phi = param$Phi, Lambda_s = param$Lambda_s, psi_s=param$psi_s,  grad = gout, hessian = hout)
  return(out)
}
