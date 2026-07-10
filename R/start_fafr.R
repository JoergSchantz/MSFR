
#' Provides some starting values for the parameters of a FR model
#'
#' This is a supporting function for \code{ecm_msfa}. The method employed is documented in the reference.
#'
#' The upper-triangular zero constraint is adopted to achieve identification,
#' as detailed in the reference, though the function can also be run without such constraint.
#' @param X_s List of lenght \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all the studies.
#' No standardization is carried out by the function.
#' @param k Number of common factors.
#' @param constraint  Constraint for ensuring identifiability. The default is "block_lower2", which
#' corresponds to the main proposal of De Vito et al. (2018). An alternative identification
#' strategy is triggered by  "block_lower1"; this is more restrictive but may work also with smaller
#' number of variables.
#' @param method Which method should be used to find the starting values? The two possibilities are \code{"adhoc"} for
#' the method described in De Vito et al. (2016), and \code{"fa"} for averaging over separate study-specific FA models.
#' Default is \code{"adhoc"}.
#' @param robust If \code{TRUE}, robust covariance matrix is used in place of the sample covariance. Default
#' is \code{FALSE}.
#' @param corr If \code{TRUE}, the analysis will employ the correlation matrix instead of the covariance matrix.
#' @param mcd If \code{TRUE}, the robust estimator used for the covariance is the same proposed in Pison et al. (2003),
#' otherwise the default value of the function \code{CovRob} of the \code{robust} library is employed. Default is
#' \code{FALSE}.
#' @return A list  containing  \code{Phi},\code{Lambda_s} and  \code{psi_s}, starting values for the model matrices.
#' @import psych
#' @export
#' @references De Vito, R., Bellio, R., Parmigiani, G. and Trippa, L. (2019). Multi-study Factor Analysis. Biometrics,  75, 337-346.
start_fr <- function(X_s, B_s, p_b, k, constraint = "block_lower2", method = "adhoc")
{
  S <- length(X_s)
  X <- Reduce('rbind', X_s)  
  B <- Reduce('rbind', B_s)
  #print(dim(B))
  #print(dim(X))
  fm1 <- lm(X ~ 0+B)
  beta = t(fm1$coefficients)
  
  X_tilde <- list()
  for(s in 1:S) X_tilde[[s]] <- X_s[[s]] - B_s[[s]]%*%t(beta)
  
  #X_used_s <- X_tilde
  X_used_s <- list()
  
  #if(corr & !robust)
  #for(s in 1:S)  X_used_s[[s]] <- scale(X_s[[s]])
  for(s in 1:S)  X_used_s[[s]] <- scale(X_tilde[[s]])
  #if(robust & corr & method=="adhoc"){
  #for(s in 1:S){
  #ogg_s <- if(mcd) covRob(X_s[[s]], estim = "mcd", quan = .75, ntrial = 1000) else covRob(X_s[[s]])
  #}
  #X_used_s[[s]] <- scale(X_s[[s]], center = ogg_s$center, scale = sqrt(diag(ogg_s$cov)))
  #}
  
  
  p <- dim(X_s[[1]])[2]
  Phi <- matrix(0, nrow=p, ncol=k)
  psi_s <- list()
  if(method=="adhoc"){
    X <- Reduce(rbind, X_used_s)
    X.pcr <- prcomp(X)
    Phi <- matrix(X.pcr$rotation[,1:k], nrow=p, ncol=k, byrow=FALSE)
    
    if (constraint == "block_lower1") {
      Phi[upper.tri(Phi)] <- 0
      for(s in 1:S){
        iniLS <- array(prcomp(X_used_s[[s]])$rotation, dim=c(p, j_s[s]))
        iniTot <- cbind(Phi, iniLS)
        iniTot[upper.tri(iniTot)] <- 0
        #Lambda_s[[s]] <-  matrix(iniTot[,(k+1):(k+j_s[s])], p , j_s[s])
        psi_s[[s]] <- fa(X_used_s[[s]], nfactors = k)$uniq
      }
    }
    
    if (constraint == "block_lower2") {
      Phi[upper.tri(Phi)] <- 0
      for(s in 1:S){
        #Lambda_s[[s]] = array(prcomp(X_used_s[[s]])$rotation, dim=c(p, j_s[s]))
        #Lambda_s[[s]][upper.tri(Lambda_s[[s]])] <- 0
        psi_s[[s]] <- fa(X_used_s[[s]], nfactors = k)$uniq
      }
    }
    
    if (constraint == "null") {
      Phi <- Phi
      for(s in 1:S){
        #Lambda_s[[s]] = array(prcomp(X_used_s[[s]])$rotation, dim=c(p, j_s[s]))
        psi_s[[s]] <- fa(X_used_s[[s]], nfactors = k)$uniq
      }
    }
  }
  #### it is important to post-process the output for avoiding sign changes
  if(method=="fa"){
    #est <- ecm_fa(X_s, tot_s = k + j_s, robust = robust, mcd = mcd, corr = corr, tol = 10^-5, nIt = 5000, trace = FALSE)
    est <- ecm_fa(X_tilde, tot_s = k, robust = robust, mcd = mcd, corr = corr, tol = 10^-5, nIt = 5000, trace = FALSE)
    Phi <- est$Omega_s[[1]][,1:k] / S
    #Lambda_s[[1]] <-  est$Omega_s[[1]][,(k+1):(k+j_s[1])]
    psi_s[[1]] <- est$psi_s[[1]]
    for(s in 2:S){
      Phi <- Phi + est$Omega_s[[s]][,1:k] / S * sign(Phi) * sign(est$Omega_s[[s]][,1:k]) ###to avoid sign changes
      #Lambda_s[[s]] <-  est$Omega_s[[s]][,(k+1):(k+j_s[s])]
      psi_s[[s]] <- est$psi_s[[s]]
    }
  }
  out <- list(Phi=Phi, psi_s=psi_s, beta=beta)
  return(out)
}
