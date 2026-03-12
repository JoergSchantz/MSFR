# TESTING FOR START_MSFA()

start_msfa <- function( X_s, B_s, p_b, k, j_s, constraint = "block_lower2", method = "adhoc" )
{
  S <- length( X_s )
  # X <- Reduce( 'rbind', X_s )
  # B <- Reduce( 'rbind', B_s )
  X <- dplyr::bind_rows( X_s )
  B <- dplyr::bind_rows( B_s )

  m1 <- stats::lm( X ~ 0 + B )
  beta <- t( fm1$coefficients )
  # X_tilde <- list()
  # for( s in 1:S ) X_tilde[[s]] <- X_s[[s]] - B_s[[s]]%*%t( beta )
  X_tilde <- lapply( 1:S, function( s ) X_s[[s]] - tcrossprod( B_s[[s]], beta ) )
  
  # X_used_s <- list()
  # for( s in 1:S ) X_used_s[[s]] <- scale( X_tilde[[s]] )
  X_used_s <- lapply( X_tilde, scale )
  
  p <- dim( X_s[[1]] )[2]
  Phi <- matrix( 0, nrow = p, ncol = k )
  Lambda_s <- psi_s <- list()
  if( method=="adhoc" ){
    # X <- Reduce( rbind, X_used_s )
    X <- dplyr::bind_rows( X_used_s )
    X.pcr <- prcomp( X )
    Phi <- matrix( X.pcr$rotation[,1:k], nrow = p, ncol = k, byrow = FALSE )
    
    if ( constraint == "block_lower1" ) {
      Phi[upper.tri( Phi )] <- 0
      for( s in 1:S ){
        iniLS <- array( prcomp( X_used_s[[s]] )$rotation, dim=c( p, j_s[s] ) )
        iniTot <- cbind( Phi, iniLS )
        iniTot[upper.tri( iniTot )] <- 0
        Lambda_s[[s]] <-  matrix( iniTot[,( k+1 ):( k+j_s[s] )], p , j_s[s] )
        psi_s[[s]] <- fa( X_used_s[[s]], nfactors = k+j_s[s] )$uniq
      }
    }
    
    if ( constraint == "block_lower2" ) {
      Phi[upper.tri( Phi )] <- 0
      for( s in 1:S ){
        Lambda_s[[s]] = array( prcomp( X_used_s[[s]] )$rotation, dim=c( p, j_s[s] ) )
        Lambda_s[[s]][upper.tri( Lambda_s[[s]] )] <- 0
        psi_s[[s]] <- fa( X_used_s[[s]], nfactors = k+j_s[s] )$uniq
      }
    }
    
    if ( constraint == "null" ) {
      Phi <- Phi
      for( s in 1:S ){
        Lambda_s[[s]] = array( prcomp( X_used_s[[s]] )$rotation, dim=c( p, j_s[s] ) )
        psi_s[[s]] <- fa( X_used_s[[s]], nfactors = k+j_s[s] )$uniq
      }
    }
  }
  #### it is important to post-process the output for avoiding sign changes
  if( method=="fa" ){
    est <- ecm_fa( X_tilde, tot_s = k + j_s, robust = robust, mcd = mcd, corr = corr, tol = 10^-5, nIt = 5000, trace = FALSE )
    Phi <- est$Omega_s[[1]][,1:k] / S
    Lambda_s[[1]] <-  est$Omega_s[[1]][,( k+1 ):( k+j_s[1] )]
    psi_s[[1]] <- est$psi_s[[1]]
    for( s in 2:S ){
      Phi <- Phi + est$Omega_s[[s]][,1:k] / S * sign( Phi ) * sign( est$Omega_s[[s]][,1:k] ) ###to avoid sign changes
      Lambda_s[[s]] <-  est$Omega_s[[s]][,( k+1 ):( k+j_s[s] )]
      psi_s[[s]] <- est$psi_s[[s]]
    }
  }
  out <- list( Phi=Phi, Lambda_s=Lambda_s, psi_s=psi_s, beta=beta )
  return( out )
}
