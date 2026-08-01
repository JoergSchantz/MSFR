
#' Provides some starting values for the parameters of a FR model
#'
#' This is a supporting function for \code{ecm_fr}. The method employed is documented in the references.
#'
#' The upper-triangular zero constraint is adopted to achieve identification,
#' as detailed in the reference, though the function can also be run without such constraint.
#' @param X_s List of lenght \eqn{S}{S}, corresponding to the number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}for all the studies.
#' No standardization is carried out by the function.
#' @param B_s List of length \eqn{S}{S}, corresponding to the number of different studies considered. 
#' Each element of the list contains a data matrix with known covariates or batch size idicators etc. 
#' @param k Number of common factors.
#' @param constraint  Constraint for ensuring identifiability. The default is "block_lower2", which
#' corresponds to the main proposal of De Vito et al. (2018). An alternative identification
#' strategy is triggered by  "block_lower1"; this is more restrictive but may work also with smaller
#' number of variables.
#' @param method Which method should be used to find the starting values? The two possibilities are \code{"adhoc"} for
#' the method described in De Vito et al. (2016), and \code{"fa"} for averaging over separate study-specific FA models.
#' Default is \code{"adhoc"}.
#' @return A list  containing  \code{Phi},\code{psi_s} and  \code{beta}, starting values for the model matrices.
#' @import psych
#' @export
#' @references De Vito, R., Bellio, R., Parmigiani, G. and Trippa, L. (2019). Multi-study Factor Analysis. Biometrics,  75, 337-346.
start_fr <- function(X_s, B_s, k, constraint = "block_lower2", method = "adhoc")
{
  S <- length(X_s)
  X <- Reduce('rbind', X_s)  
  B <- Reduce('rbind', B_s)
  fm1 <- stats::lm(X ~ 0+B)
  beta = t(fm1$coefficients)
  
  X_tilde <- vector("list", S)
  for(s in 1:S) X_tilde[[s]] <- X_s[[s]] - tcrossprod(B_s[[s]], beta)
  
  X_used_s <- vector("list", S)
  for(s in 1:S)  X_used_s[[s]] <- scale(X_tilde[[s]])
  
  p <- dim(X_s[[1]])[2]
  Phi <- matrix(0, nrow=p, ncol=k)
  psi_s <- vector("list", S)
  if(method=="adhoc"){
    X <- Reduce(rbind, X_used_s)
    X.pcr <- stats::prcomp(X)
    Phi <- matrix(X.pcr$rotation[,1:k], nrow=p, ncol=k, byrow=FALSE)
    
    if (constraint == "block_lower1") {
      Phi[upper.tri(Phi)] <- 0
      for(s in 1:S){
        iniLS <- array(stats::prcomp(X_used_s[[s]])$rotation, dim=c(p, j_s[s]))
        iniTot <- cbind(Phi, iniLS)
        iniTot[upper.tri(iniTot)] <- 0
        psi_s[[s]] <- psych::fa(X_used_s[[s]], nfactors = k)$uniq
      }
    }
    
    if (constraint == "block_lower2") {
      Phi[upper.tri(Phi)] <- 0
      for(s in 1:S){
        psi_s[[s]] <- psych::fa(X_used_s[[s]], nfactors = k)$uniq
      }
    }
    
    if (constraint == "null") {
      for(s in 1:S){
        psi_s[[s]] <- psych::fa(X_used_s[[s]], nfactors = k)$uniq
      }
    }
  }
  #### it is important to post-process the output for avoiding sign changes
  if(method=="fa"){
    est <- ecm_fa(X_tilde, tot_s = k, tol = 10^-5, nIt = 5000, trace = FALSE)
    Phi <- est$Omega_s[[1]][,1:k] / S
    psi_s[[1]] <- est$psi_s[[1]]
    for(s in 2:S){
      Phi <- Phi + est$Omega_s[[s]][,1:k] / S * sign(Phi) * sign(est$Omega_s[[s]][,1:k]) ###to avoid sign changes
      psi_s[[s]] <- est$psi_s[[s]]
    }
  }
  out <- list(Phi=Phi, psi_s=psi_s, beta=beta)
  return(out)
}
