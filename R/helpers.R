# This is Joerg's playground for the MSFR package 
## Probably will result in a helper file, which contains all 'internal' functionalities
## lets goooo

#' @keywords internal
#' Computes trace of a matrix
.tr <- function( A ) sum( diag( A ) )

#' @keywords internal
#' Inverts often used Psi matrix more efficiently than solve()
.inv_Psi <- function( Psi ) diag( 1 / diag( Psi ) )

#' @keywords internal
#' Same as .inv_Psi(), but returns the length-p diagonal as a plain vector
#' instead of a dense p x p matrix.
.inv_Psi_vec <- function( Psi ) 1 / diag( Psi )



#' @keywords internal
.loglik_ecm <- function( Sig_s1, ds_s, n_s, cov_s )
{
  S <- length( n_s )
  #####log likelihood value for each study
  val_s <- sapply(
    1:S,
    function( s ) {
      - ( n_s[s]/2 ) * log( ds_s[[s]] ) - ( n_s[s]/2 ) * .tr( Sig_s1[[s]] %*% cov_s[[s]] )
    }
  )
  #####sum of each study-likelihood
  val_tot <- sum( val_s )
  return( val_tot )
}


#' @keywords internal
.param2vect <- function(param, constraint)
{
  p <- length(param$Psi_s[[1]])
  S <- length(param$Psi_s)
  phi_vals <- as.vector(param$Phi[lower.tri(param$Phi, diag = TRUE)])
  k <- ncol(param$Phi)
  lambda_vals <- psi_vals <- c()
  # j_s <- c()
  
  if(constraint=="block_lower1"){
    for(s in 1:S){
      Lambda_s <- param$Lambda_s[[s]][-(1:k),]
      lambda_vals <- c(lambda_vals, as.vector(Lambda_s[lower.tri(Lambda_s, diag = TRUE)]))
      psi_vals <- c(psi_vals, param$Psi_s[[s]])
      #j_s[[s]] <- ncol(param$Lambda_s[[s]])
    }
  }
  if(constraint=="block_lower2"){
    for(s in 1:S){
      Lambda_s <- param$Lambda_s[[s]]
      lambda_vals <- c(lambda_vals, as.vector(Lambda_s[lower.tri(Lambda_s, diag = TRUE)]))
      psi_vals <- c(psi_vals, param$Psi_s[[s]])
      #j_s[[s]] <- ncol(param$Lambda_s[[s]])
    }
  }
  theta <- c(phi_vals, lambda_vals, psi_vals)
  return(theta)
}


#' @keywords internal
.vect2param <- function(vect, param.struct, constraint, p, k, j_s)
{
  S <- length(j_s)
  nP <- k * p - k * ( k - 1) / 2
  if (constraint == "block_lower1")  { nL <- j_s * (p - k)  - j_s *  (j_s - 1) / 2}
  if (constraint == "block_lower2")  { nL <- j_s * p  - j_s *  (j_s - 1) / 2}
  phi_vals <- vect[1:nP]
  Phi <- matrix(0, nrow=p, ncol=k)
  Phi[lower.tri(Phi, diag = TRUE)] <- phi_vals
  Lambda_s <- param.struct$Lambda_s
  psi_s <- param.struct$psi_s
  for(s in 1:S){
    nL_s  <- if(s==1) 0 else sum(nL[1:(s-1)])
    ind <-  (nP + nL_s + 1):(nP + nL_s + nL[s])
    lambda_vals_s <-  vect[ind]
    Lambda_s[[s]][Lambda_s[[s]]!=0] <- lambda_vals_s
    ind_s <- (nP + sum(nL) + p * (s-1) + 1):(nP + sum(nL) + p * s)
    psi_vals_s <- vect[ind_s]
    psi_s[[s]] <-  psi_vals_s
  }
  return(list(Phi=Phi, Lambda_s=Lambda_s, psi_s=psi_s))
}


#' @keywords internal
#### loglikelihood function re-expressed as a function of the model parameters
#### theta: c(Phi, Lambda_1,..,Lambda_S,Psi_1,..,Psi_S)
.loglik_int <- function(theta, n_s, cov_s, k, j_s, constraint)
{
  S <- length(n_s)
  p <- ncol(cov_s[[1]])
  nP <- k * p - k * ( k - 1) / 2
  if (constraint == "block_lower1")  { nL <- j_s * (p - k)  - j_s *  (j_s - 1) / 2}
  if (constraint == "block_lower2")  { nL <- j_s * p  - j_s *  (j_s - 1) / 2}
  phi_vals <- theta[1:nP]
  Phi <- matrix(0, p, k)
  Phi[lower.tri(Phi, diag = TRUE)] <- phi_vals
  out <- 0
  for(s in 1:S){
    nL_s  <- if(s==1) 0 else sum(nL[1:(s-1)])
    ind <-  (nP + nL_s + 1):(nP + nL_s + nL[s])
    omega_vals_s <- c(phi_vals, theta[ind])
    Omega_s <- matrix(0, p, k + j_s[s])
    if (constraint == "block_lower1")  Omega_s[lower.tri(Omega_s, diag = TRUE)] <- omega_vals_s
    if (constraint == "block_lower2")
    {
      Lambda_s <- matrix(0, p, j_s[s])
      Lambda_s[lower.tri(Lambda_s, diag = TRUE)] <- theta[ind]
      Omega_s <- cbind(Phi, Lambda_s)
    }
    ind_s <- (nP + sum(nL) + p * (s-1) + 1):(nP + sum(nL) + p * s)
    psi_vals_s <- theta[ind_s]
    Psi_s1 <- diag(1 / psi_vals_s)
    D1L_s <- statmod::vecmat(1 / sqrt(psi_vals_s), Omega_s)
    LDL_s <- crossprod(D1L_s)
    A <- diag(k + j_s[s]) + LDL_s
    A1 <- chol2inv(chol(A))
    D2L_s <- statmod::vecmat(1/psi_vals_s, Omega_s)
    Sig_s1 <- Psi_s1 - D2L_s  %*% A1 %*% t(D2L_s)
    log_ds_s <-  log(det(A)) + sum(log(psi_vals_s))
    out  <- out - (n_s[s]/2) * log_ds_s - (n_s[s]/2) * .tr(Sig_s1 %*% cov_s[[s]])
  }
  return(out)
}