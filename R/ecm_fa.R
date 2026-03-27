#' Estimates the parameters of study-specific FA models
#'
#' Maximum likelihood estimation of study-specific FA models parameters via the ECM
#' algorithm, adopting the upper-triangular zero constraint to achieve identification
#' for each loading matrix. Note: the function can also estimate a FA model for a single
#' study, by specifiyng \code{X_s = list(data)}, where \code{data} is the data matrix.
#' @param X_s List of lenght \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all the studies.
#' @param tot_s Number of latent factors for each study. A vector of positive integers of length \eqn{S}{S}.
#' @param nIt Maximum number of iterations for the ECM algorithm. Default is 50000.
#' @param tol Tolerance for declaring convergence of the ECM algorithm. Default is 10^-7.
#' @param block_lower Should the upper-triangular zero constraint be enforced? Default is \code{TRUE}
#' (strongly suggested).
#' @param robust If \code{TRUE}, robust covariance matrix is used in place of the sample covariance. Default
#' is \code{FALSE}.
#' @param corr If \code{TRUE}, the analysis will employ the correlation matrix instead of the covariance matrix.
#' @param mcd If \code{TRUE}, the robust estimator used for the covariance is the same proposed in Pison et al. (2003),
#' otherwise the default value of the function \code{CovRob} of the \code{robust} library is employed. Default is
#' \code{FALSE}.
#' @param trace If \code{TRUE} then trace information is being printed every \code{traceIT} iterations of the ECM algorithm.
#' @param traceIT Frequency of tracing information.
#' @return A list  containing the following components:
#' \item{\code{Omega_s}, \code{Psi_s}}{the estimated model matrices.}
#' \item{loglik}{the value of the log likelihood function at the final estimates.}
#' \item{\code{AIC, BIC}}{model selection criteria at the estimate.}
#' \item{\code{npar}}{number of model parameters.}
#' \item{iter}{the number of ECM iterations performed.}
#' @export
#' @import robust psych
#' @references De Vito, R., Bellio, R., Trippa, L. and Parmigiani, G. (2019). Multi-study Factor Analysis. Biometrics,  75, 337-346.
#' @references Pison, G., Rousseeuw, P.J., Filzmoser, P. and Croux, C. (2003). Robust factor analysis. Journal
#' Multivariate Analysis, 84, 145-172.
ecm_fa <- function(X_s, tot_s, nIt = 50000, tol = 10^-7, block_lower = TRUE, robust = FALSE, corr = TRUE, mcd = FALSE, trace = TRUE, traceIT = 1000)
{
  Omega_s <- list()
  Psi_s <- psi_s <- list()
  #######
  p <- ncol(X_s[[1]])
  S <- length(X_s)
  n_s <- numeric(S)
  #######defining objects
  Psi_s1 <- list()
  cov_s <- list()
  Phi <- matrix(0, nrow=p, ncol=1)
  ######1st round of cycle
  for(s in 1:S){
    n_s[s] <-  dim(X_s[[s]])[[1]]
    if((!robust) & (!corr)) cov_s[[s]] <- cov(X_s[[s]])
    if((!robust) & corr) cov_s[[s]] <- cor(X_s[[s]])
    if(robust & mcd) cov_s[[s]] <- covRob(X_s[[s]], estim = "mcd", quan = .75, ntrial = 1000, corr = corr)$cov
    if(robust & (!mcd)) cov_s[[s]] <- covRob(X_s[[s]], corr = corr)$cov
    FA.s <- factanal(X_s[[s]], factors = tot_s[[s]], covmat = cov_s[[s]],  n.obs=nrow(X_s[[s]]), rotation = "none")
    Omega_s[[s]] <- FA.s$loadings
    Psi_s[[s]] <- diag(FA.s$uniq)
    Psi_s1[[s]] <-  diag(1/diag(Psi_s[[s]]))
  }
  ######E-step
  out <- .exp_values(Phi, Omega_s, Psi_s, Psi_s1, cov_s, getdet = TRUE)
  Sig_s1 <- out$Sig_s1; ds_s= out$ds_s;
  l_stop0 <- 0
  lm1 <- 0
  l0 <- .loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
  for (i in (1:nIt))
  {
    ###########CM1 ---------------------------------------------------------------------------------------

    ######expected values
    out <- .exp_values(Phi, Omega_s, Psi_s, Psi_s1, cov_s)
    Txsfs <- out$Txsfs; Txsfcs <- out$Txsfcs; Tfsfs <- out$Tfsfs; Tfcsfcs <- out$Tfcsfcs; Tfcsfs <- out$Tfcsfs
    ######update  of Phi_s
    Psi_new <- list()
    Psi_new1 <- list()

    for(s in 1:S){
      Psi_new[[s]]  <- diag(cov_s[[s]] + Omega_s[[s]] %*% Tfsfs[[s]] %*% t(Omega_s[[s]]) -  2*Txsfs[[s]] %*% t(Omega_s[[s]]) )
      Psi_new[[s]] <- diag(Psi_new[[s]])
      ##########inverse
      Psi_new1[[s]] <- diag(1/diag(Psi_new[[s]]))
    }

    ###########CM2 ---------------------------------------------------------------------------------------

    ######expected values
    out<- .exp_values(Phi, Omega_s, Psi_new, Psi_new1, cov_s)
    Txsfs <- out$Txsfs; Txsfcs <- out$Txsfcs; Tfsfs <- out$Tfsfs;
    Tfcsfcs <- out$Tfcsfcs; Tfcsfs <- out$Tfcsfs

    ######update of Phi: not needed


    ########CM3 ---------------------------------------------------------------------------------------

    ######expected values
    out <- .exp_values(Phi, Omega_s, Psi_new, Psi_new1, cov_s)
    Txsfs <- out$Txsfs; Txsfcs <- out$Txsfcs; Tfsfs <- out$Tfsfs;
    Tfcsfcs <-  out$Tfcsfcs; Tfcsfs <- out$Tfcsfs

    ######update of Phi
    Omega_new <- list()
    for(s in 1:S) {
      Omega_new[[s]] <- matrix((Txsfs[[s]] %*% solve(Tfsfs[[s]])), p, tot_s[s])
      Omega_new[[s]][upper.tri(Omega_new[[s]])] <- 0
    }

    ###########stopping rule
    out <- .exp_values(Phi, Omega_new, Psi_new, Psi_new1, cov_s, getdet = TRUE)

     Sig_s1 <- out$Sig_s1
    ds_s <- out$ds_s
    l1 <- .loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
    a <- (l1 - l0)/ (l0-lm1)
    l_stop <- lm1 + (1/ (1-a)) * (l0-lm1)
    l0 <- .loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
    if((trace) & (i %% 100 == 0))  cat("i=", i, "Criterion for convergence ", abs(l_stop-l_stop0), "\n")
    if( (abs(l_stop-l_stop0)<tol) & i > 1 & l_stop != Inf) break
    Psi_s <- Psi_new
    Omega_s <- Omega_new
    Psi_s1 <- Psi_new1
    lm1 <- l0
    l0 <- l1
    l_stop0 <- l_stop
  }
  ############return output
  for(s in 1:S) psi_s[[s]] <- diag(Psi_s[[s]])
  npar <- p * S + sum(tot_s * (p - (tot_s - 1) / 2))
  n_tot <- sum(n_s)
  AIC <- -2 * l1 + npar * 2
  BIC <- -2 * l1 + npar * log(n_tot)
  res <- list(Omega_s = Omega_s, psi_s = psi_s, loglik = l1, AIC = AIC,
              BIC = BIC, npar=npar, iter = i)
  return(res)
}
