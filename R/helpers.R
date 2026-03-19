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
#' Implementation of Woodbury Identity. This will most likely only work 
#' inside this package and is not meant to be used outside of it.
#' @references O. Morgenstern and M. A. Woodbury, “Stability of Inverses of Input-Output Matrices,” Econometrica 18 (1950): 190.
.wb_identity <- function( A, W, I) {
  # A: being either Lambda_s[[s]] or Phi matrix
  # I: Study (un-)specific Identity matrix
  cp <- crossprod( A, W )               # to save computation time
  solve( I + cp %*% A ) %*% cp
}


#' @importFrom statmod vecmat 
#' @importFrom statmod matvec
#' @keywords internal
.exp_values <- function(Phi, Lambda_s, Psi_s, Psi_s1, cov_s, X_s_tilde, getdet = FALSE)
{
  k <- dim(Phi)[2]
  I_k <- diag(1, k)
  S <- length(Lambda_s)
  
  ###defining objects
  j_s <- numeric(S)
  I_j <- list()
  Sig_s <- list()
  ds_s <- list()
  I_tot <- list()
  LambTOT <- list()
  Sig_s1 <- list()
  delta_Lambda <- list()
  delta_Phi <- list()
  Delta_Lambda <- list()
  Delta_Phi <- list()
  Covfcfs <- list()
  Txsfs <- list()
  Txsfcs <- list()
  Tfsfs <- list()
  Tfcsfcs <- list()
  Tfcsfs <- list()
  wb1_f <- list()
  wb1_l <- list()
  wb2_f <- list()
  wb2_l <- list()
  E_fis_x_is <- list()
  E_lis_x_is <- list()
  
  for (s in 1:S){
    ds_s[[s]] <- NULL
    j_s[s] <- c( dim( Lambda_s[[s]] )[[2]] )
    I_j[[s]] <- diag( 1, j_s[s] )
    Sig_s[[s]] <- tcrossprod( Phi ) + tcrossprod( Lambda_s[[s]] ) + Psi_s[[s]]
    if (getdet) { ds_s[[s]] <- det( Sig_s[[s]] ) }
    I_tot[[s]] <- diag( 1, k + j_s[s] )
    LambTOT[[s]] <- cbind( Phi, Lambda_s[[s]] )
    # needs more work here ... seems questionable in general (see explanation in paper)
    Sig_s1[[s]] <- Psi_s1[[s]] - (statmod::vecmat(diag(Psi_s1[[s]]), LambTOT[[s]]) %*%
                                    solve(I_tot[[s]] + (t(LambTOT[[s]]) %*% statmod::vecmat(diag(Psi_s1[[s]]),
                                                                                            LambTOT[[s]]))) %*% statmod::matvec(t(LambTOT[[s]]), diag(Psi_s1[[s]])))
    # ----
    ##new part for beta
    inv_Psi_s <- .inv_Psi( Psi_s[[s]] ) 
    wb1_f[[s]] <- inv_Psi_s - inv_Psi_s %*% Lambda_s[[s]] %*% .wb_identity( Lambda_s[[s]], inv_Psi_s, I_j[[s]] )
    wb1_l[[s]] <- inv_Psi_s - inv_Psi_s %*% Phi %*% .wb_identity( Phi, inv_Psi_s, I_k ) 
    wb2_f[[s]] <- .wb_identity( Phi, wb1_f[[s]], I_k ) 
    wb2_l[[s]] <- .wb_identity( Lambda_s[[s]], wb1_l[[s]], I_j[[s]] ) 
    E_fis_x_is[[s]] <- tcrossprod( wb2_f[[s]], X_s_tilde[[s]] )
    E_lis_x_is[[s]] <- tcrossprod( wb2_l[[s]], X_s_tilde[[s]] )

    # delta_Lambda[[s]] <- crossprod( Lambda_s[[s]], Sig_s1[[s]] )
    # delta_Phi[[s]] <- crossprod( Phi, Sig_s1[[st]] )

    # THE FOLLOWING SUBSTITUTIONS HAVE NOT BEEN TESTED YET BUT SEEM REASONABLE ----
    # Delta_Lambda[[s]] <- I_j[[s]] - (t(Lambda_s[[s]]) %*% Sig_s1[[s]] %*% Lambda_s[[s]])
    Delta_Lambda[[s]] <- I_j[[s]] - ( wb2_l[[s]] %*% Lambda_s[[s]] )
    # Delta_Phi[[s]] <- I_k - (t(Phi) %*% Sig_s1[[s]] %*% Phi)
    Delta_Phi[[s]] <- I_k - ( wb2_f[[s]] %*% Phi )
    # Covfcfs[[s]] <- -t(Phi) %*% Sig_s1[[s]] %*% Lambda_s[[s]]
    Covfcfs[[s]] <- -wb2_f[[s]] %*% Lambda_s[[s]]
    # Txsfs[[s]] <- cov_s[[s]] %*% t(delta_Lambda[[s]])
    Txsfs[[s]] <- tcrossprod( cov_s[[s]], wb2_l[[s]] )
    # Txsfcs[[s]] <- cov_s[[s]] %*% t(delta_Phi[[s]])
    Txsfcs[[s]] <- tcrossprod( cov_s[[s]], wb2_f[[s]] )
    # Tfsfs[[s]] <- delta_Lambda[[s]] %*% cov_s[[s]] %*% t(delta_Lambda[[s]]) + Delta_Lambda[[s]]
    Tfsfs[[s]] <- wb2_l[[s]] %*% tcrossprod( cov_s[[s]], wb2_l[[s]] ) + Delta_Lambda[[s]]
    # Tfcsfcs[[s]] <- delta_Phi[[s]] %*% cov_s[[s]] %*% t(delta_Phi[[s]]) + Delta_Phi[[s]]
    Tfcsfcs[[s]] <- wb2_f[[s]] %*% tcrossprod( cov_s[[s]], wb2_f[[s]] ) + Delta_Phi[[s]]
    # Tfcsfs[[s]] <- delta_Phi[[s]] %*% cov_s[[s]] %*% t(delta_Lambda[[s]]) + Covfcfs[[s]]
    Tfcsfs[[s]] <- wb2_f[[s]] %*% tcrossprod( cov_s[[s]], wb2_l[[s]] ) + Covfcfs[[s]]
    # ----
  }
  return(list(Txsfs = Txsfs, Txsfcs = Txsfcs, Tfsfs = Tfsfs,
              Tfcsfcs =  Tfcsfcs, Tfcsfs = Tfcsfs, 
              E_fis_x_is = E_fis_x_is, E_lis_x_is = E_lis_x_is, 
              ds_s=ds_s,  Sig_s1 = Sig_s1))
}


#' @importFrom statmod vecmat 
#' @importFrom statmod matvec
#' @keywords internal
.exp_values_fr <- function( Phi, Psi_s, Psi_s1, cov_s, X_s_tilde, getdet = FALSE )
{
  k <- dim( Phi )[2]
  I_k <- diag( 1, k )
  S <- length( Psi_s )
  
  ###defining objects
  Sig_s <- list()
  ds_s <- list()
  I_tot <- list()
  # LambTOT <- list()
  Sig_s1 <- list()
  delta_Phi <- list()
  Delta_Phi <- list()
  Covfcfs <- list()
  Txsfcs <- list()
  Tfcsfcs <- list()
  Woodbury_f <- list()
  E_fis_x_is <- list()
  
  for ( s in 1:S ){
    ds_s[[s]] <- NULL
    #j_s[s] <- c(dim(Lambda_s[[s]])[[2]])
    # Sig_s[[s]] <- Phi %*% t(Phi) + Psi_s[[s]]
    Sig_s[[s]] <- tcrossprod( Phi ) + Psi_s[[s]]
    if ( getdet ) ds_s[[s]] <- det( Sig_s[[s]] )
    # I_tot[[s]] <- diag( 1, k )
    # LambTOT[[s]] <- Phi
    # Sig_s1[[s]] <- Psi_s1[[s]] - (statmod::vecmat(diag(Psi_s1[[s]]), LambTOT[[s]]) %*%
    #                                solve(I_tot[[s]] + (t(LambTOT[[s]]) %*% statmod::vecmat(diag(Psi_s1[[s]]),
    #
    Sig_s1[[s]] <- Psi_s1[[s]] - (statmod::vecmat(diag(Psi_s1[[s]]), Phi) %*%
                                    solve(I_k + (t(Phi) %*% statmod::vecmat(diag(Psi_s1[[s]]),
                                                                                             Phi ))) %*% statmod::matvec(t( Phi ), diag(Psi_s1[[s]])))
    # delta_Phi[[s]] <- t(Phi) %*% Sig_s1[[s]]
    delta_Phi[[s]] <- crossprod( Phi, Sig_s1[[s]] )
    # Delta_Phi[[s]] <- I_k - (t(Phi) %*% Sig_s1[[s]] %*% Phi)
    Delta_Phi[[s]] <- I_k - ( crossprod( Phi, Sig_s1[[s]] ) %*% Phi )
    # Txsfcs[[s]] <- cov_s[[s]] %*% t(delta_Phi[[s]])
    Txsfcs[[s]] <- tcrossprod( cov_s[[s]], delta_Phi[[s]] )
    # Tfcsfcs[[s]] <- delta_Phi[[s]] %*% cov_s[[s]] %*% t(delta_Phi[[s]]) + Delta_Phi[[s]]
    Tfcsfcs[[s]] <- delta_Phi[[s]] %*% tcrossprod( cov_s[[s]], delta_Phi[[s]] ) + Delta_Phi[[s]]

    ##new part for beta
    # Woodbury_f[[s]] <- solve(I_k + t(Phi) %*%  diag(1/diag(Psi_s[[s]])) %*% Phi) %*% t(Phi) %*% diag(1/diag(Psi_s[[s]]))
    Woodbury_f[[s]] <- .wb_identity( Phi, .inv_Psi( Psi_s[[s]] ), I_k )
    # E_fis_x_is[[s]] <- Woodbury_f[[s]] %*% t(X_s_tilde[[s]])
    E_fis_x_is[[s]] <- tcrossprod( Woodbury_f[[s]], X_s_tilde[[s]] )

  }
  return(
    list( 
      Txsfcs = Txsfcs,
      Tfcsfcs =  Tfcsfcs,
      E_fis_x_is = E_fis_x_is, 
      ds_s=ds_s,
      Sig_s1 = Sig_s1
    )
  )
}


#' @keywords internal
.loglik_ecm <- function( Sig_s1, ds_s, n_s, cov_s )
{
  S <- length( n_s )
  #####log likelihood value for each study
  # val_s <- c()
  # for(s in 1:S){
  #   val_s[s] <- - (n_s[s]/2) * log(ds_s[[s]]) - (n_s[s]/2) * .tr(Sig_s1[[s]] %*% cov_s[[s]])
  # }
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
  p <- length(param$psi_s[[1]])
  S <- length(param$psi_s)
  phi_vals <- as.vector(param$Phi[lower.tri(param$Phi, diag = TRUE)])
  k <- ncol(param$Phi)
  lambda_vals <- psi_vals <- j_s <- c()
  
  if(constraint=="block_lower1"){
    for(s in 1:S){
      Lambda_s <- param$Lambda_s[[s]][-(1:k),]
      lambda_vals <- c(lambda_vals, as.vector(Lambda_s[lower.tri(Lambda_s, diag = TRUE)]))
      psi_vals <- c(psi_vals, param$psi_s[[s]])
      j_s[[s]] <- ncol(param$Lambda_s[[s]])}
  }
  if(constraint=="block_lower2"){
    for(s in 1:S){
      Lambda_s <- param$Lambda_s[[s]]
      lambda_vals <- c(lambda_vals, as.vector(Lambda_s[lower.tri(Lambda_s, diag = TRUE)]))
      psi_vals <- c(psi_vals, param$psi_s[[s]])
      j_s[[s]] <- ncol(param$Lambda_s[[s]])}
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