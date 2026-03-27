#' @keywords internal
.get_exp_xl <- function( delta_Lambda, cov_s ) {
  tcrossprod( cov_s, delta_Lambda )
}

#' @keywords internal
.get_exp_xf <- function( delta_Phi, cov_s ) {
  tcrossprod( cov_s, delta_Phi )
}

#' @keywords internal
.get_exp_ll <- function( Lambda_s, delta_Lambda, cov_s, I_j ) {
  delta_Lambda %*% tcrossprod( cov_s, delta_Lambda ) + ( I_j - ( delta_Lambda %*% Lambda_s ) )
}

#' @keywords internal
.get_exp_ff <- function( Phi, delta_Phi, cov_s, I_k ) {
  delta_Phi %*% tcrossprod( cov_s, delta_Phi ) + ( I_k - ( delta_Phi %*% Phi ) )
}

#' @keywords internal
.get_exp_fl <- function( Lambda_s, delta_Phi, delta_Lambda, cov_s ) {
  delta_Phi %*% tcrossprod( cov_s, delta_Lambda ) + ( -delta_Phi %*% Lambda_s )
}

#' @keywords internal
.get_exp_f <- function( delta_Phi, X_s_tilde ) {
  tcrossprod( delta_Phi, X_s_tilde )
}

#' @keywords internal
.get_exp_l <- function( delta_Lambda, X_s_tilde ) {
  tcrossprod( delta_Lambda, X_s_tilde )
}

#' @keywords internal
#' perfroms calculations in CM step 1 when .exp_values() is called
.step_cm_1 <- function( Phi, Lambda_s, delta_Phi, delta_Lambda, cov_s, I_k, I_j ) {
  exp_xl <- Map( .get_exp_xl, delta_Lambda, cov_s )
  exp_xf <- Map( .get_exp_xf, delta_Phi, cov_s )
  exp_ll <- Map( .get_exp_ll, Lambda_s, delta_Lambda, cov_s, I_j )
  exp_ff <- Map( .get_exp_ff, list(Phi), delta_Phi, cov_s, list(I_k) )
  exp_fl <- Map( .get_exp_fl, Lambda_s, delta_Phi, delta_Lambda, cov_s )

  return(
    list(
      exp_xl = exp_xl, 
      exp_xf = exp_xf, 
      exp_ll = exp_ll, 
      exp_ff = exp_ff, 
      exp_fl = exp_fl
    )
  )
}

#' @keywords internal
#' same as .step_cm_1()
.step_cm_2 <- function( Phi, Lambda_s, delta_Phi, delta_Lambda, cov_s, I_k ) {
  exp_ff <- Map( .get_exp_ff, list(Phi), delta_Phi, cov_s, list(I_k) )
  exp_fl <- Map( .get_exp_fl, Lambda_s, delta_Phi, delta_Lambda, cov_s )
  exp_xf <- Map( .get_exp_xf, delta_Phi, cov_s )

  return(
    list(
      exp_xf = exp_xf,
      exp_ff = exp_ff, 
      exp_fl = exp_fl
    )
  )
}

#' @keywords internal
.step_cm_3 <- function( Lambda_s, delta_Phi, delta_Lambda, cov_s, I_j ) {
  exp_xl <- Map( .get_exp_xl, delta_Lambda, cov_s )
  exp_ll <- Map( .get_exp_ll, Lambda_s, delta_Lambda, cov_s, I_j )
  exp_fl <- Map( .get_exp_fl, Lambda_s, delta_Phi, delta_Lambda, cov_s )

  return(
    list(
      exp_xl = exp_xl, 
      exp_ll = exp_ll, 
      exp_fl = exp_fl
    )
  )
}

#' @keywords internal
.step_cm_4 <- function( delta_Phi, delta_Lambda, X_s_tilde ) {
  exp_l <- Map( .get_exp_l, delta_Lambda, X_s_tilde )
  exp_f <- Map( .get_exp_f, delta_Phi, X_s_tilde )

  return(
    list(
      exp_f = exp_f, 
      exp_l = exp_l 
    )
  )
}

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
#' @keywords internal
#' 2nd application of the Woodbury Identity
.wb_identity2 <- function( inv_Psi, A, I ){
  inv_Psi - inv_Psi %*% A %*% .wb_identity( A, inv_Psi, I )
}

#' @keywords internal
.exp_values <- function(
  Phi, 
  Lambda_s, 
  Psi_s,  
  CM_step = 1, # out of (1,2,3,4)
  cov_s = NULL, # not needed for cm_step == 4
  X_s_tilde = NULL # only needed for cm_step == 4
)
{
  k <- dim( Phi )[2]
  I_k <- diag( 1, k )
  I_j <- lapply( Lambda_s, function( s ) diag( 1, dim( s )[2] ) )

  inv_Psi_s <- lapply( Psi_s, .inv_Psi )
  
  wb_f <- Map( .wb_identity2, inv_Psi_s, Lambda_s, I_j )
  wb_l <- Map( .wb_identity2, inv_Psi_s, list(Phi), list(I_k) )
  
  delta_Phi <- Map( .wb_identity, list(Phi), wb_f, list(I_k) )
  delta_Lambda <-  Map( .wb_identity, Lambda_s, wb_l, I_j )

  if ( CM_step == 1 ) {
    return(
      .step_cm_1( Phi, Lambda_s, delta_Phi, delta_Lambda, cov_s, I_k, I_j )
    )
  }

  if ( CM_step == 2 ) {
    return(
      .step_cm_2( Phi, Lambda_s, delta_Phi, delta_Lambda, cov_s, I_k )
    )
  }

  if ( CM_step == 3 ) {
    return(
      .step_cm_3( Lambda_s, delta_Phi, delta_Lambda, cov_s, I_j )
    )
  }

  if ( CM_step == 4 ) {
    return( 
      .step_cm_4( delta_Phi, delta_Lambda, X_s_tilde )
    )
  }
}