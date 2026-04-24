#' @keywords internal
#' --- helper: sparse scaling ---
scale_sparse <- function( M ) {
  cm <- Matrix::colMeans( M )
  M_centered <- sweep( M, 2, cm, "-" )

  csd <- sqrt( Matrix::colMeans( M_centered^2 ) )
  csd[csd == 0] <- 1

  M_scaled <- sweep( M_centered, 2, csd, "/" )
  return( M_scaled )
}
#' @keywords internal
#' --- helper: sparse PCA via truncated SVD ---
pca_sparse <- function( M, k ) {
  sv <- irlba::irlba( M, nv = k )
  return( sv$v )
}
#' @keywords internal
#' --- helper: sparse covariance ---
cov_sparse <- function( M ) {
  n <- nrow( M )
  Matrix::crossprod( M ) / ( n - 1 )
}
#' @keywords internal
#' --- helper: approximate uniquenesses ---
psi_from_svd <- function( Smat, r ) {
  sv <- irlba::irlba( Smat, nv = r )
  approx <- sv$v %*% diag( sv$d ) %*% t( sv$v )
  psi <- diag( Smat - approx )
  psi[psi < 1e-8] <- 1e-8
  return( psi )
}

#'@export
start_msfa_sparse <- function( X_s, B_s, p_b, k, j_s,
                             constraint = "block_lower2",
                             method = "adhoc" ) {

  if ( !requireNamespace( "Matrix", quietly = TRUE ) )
    stop( "Package 'Matrix' required" )
  if ( !requireNamespace( "irlba", quietly = TRUE ) )
    stop( "Package 'irlba' required" )

  S <- length( X_s )



  ## --- stack data ( still sparse ) ---
  X <- do.call( rbind, X_s )
  B <- do.call( rbind, B_s )

  ## --- sparse regression: replaces lm(  ) ---
  BtB <- crossprod(B)
  BtX <- crossprod(B, X)
  beta <- solve(BtB, BtX)

  ## --- residuals ---
  X_tilde <- Map(
    function( X, B, beta ){
      X - B %*% beta
    },
    X_s, B_s, list( beta )
  )


  ## --- scaling ---
  X_used_s <- lapply( X_tilde, scale_sparse )

  p <- ncol( X_s[[1]] )
  Phi <- Matrix::Matrix( 0, nrow = p, ncol = k, sparse = FALSE )
  Lambda_s <- vector( "list", S )
  psi_s <- vector( "list", S )

  if ( method == "adhoc" ) {

    ## --- global PCA ---
    X_all <- do.call( rbind, X_used_s )
    Phi <- pca_sparse( X_all, k )

    if ( constraint %in% c( "block_lower1", "block_lower2" ) ) {
      Phi[upper.tri( Phi )] <- 0
    }

    for ( s in 1:S ) {

      ## --- study-specific PCA ---
      Ls <- pca_sparse( X_used_s[[s]], j_s[s] )

      if ( constraint == "block_lower1" ) {
        iniTot <- cbind( Phi, Ls )
        iniTot[upper.tri( iniTot )] <- 0
        Lambda_s[[s]] <- iniTot[, ( k + 1 ):( k + j_s[s] ), drop = FALSE]
      }

      if ( constraint == "block_lower2" ) {
        Ls[upper.tri( Ls )] <- 0
        Lambda_s[[s]] <- Ls
      }

      if ( constraint == "null" ) {
        Lambda_s[[s]] <- Ls
      }

      ## --- approximate uniqueness ---
      S_s <- cov_sparse( X_used_s[[s]] )
      psi_s[[s]] <- psi_from_svd( S_s, k + j_s[s] )
    }
  }

  if ( method == "fa" ) {
    stop( "method='fa' not supported in sparse version ( requires dense FA )" )
  }

  return( list( 
    Phi = Phi,
    Lambda_s = Lambda_s,
    psi_s = psi_s,
    beta = t(beta)
   ) )
}