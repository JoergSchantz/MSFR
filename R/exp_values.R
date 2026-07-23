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
#'
#' \code{W} is the "D^-1" term of the Woodbury identity, i.e. it's always
#' symmetric p x p. In MSFR's double reduction (see \code{.exp_values()}),
#' the *second* application's \code{W} is itself an already Woodbury-reduced
#' block (\code{wb_f}/\code{wb_l}) and is genuinely dense. But in every
#' *first* application -- which is the only one FA/FR ever need, since
#' neither has a second factor block to reduce against -- \code{W} is always
#' \code{Psi_s^-1}, which is diagonal. Passing that as a length-p vector
#' (e.g. via \code{.inv_Psi_vec()}) instead of a dense p x p matrix (via
#' \code{.inv_Psi()}) lets both functions skip the O(p^2) allocation and
#' replace an O(p^2 q) dense multiply with an O(p q) row-scaling one; the
#' final result is still assembled as a dense matrix where required (i.e. in
#' \code{.wb_identity2()}), since \code{Sig_s^-1} genuinely isn't diagonal.
#' Passing a dense matrix still works exactly as before -- this only adds a
#' faster path, it doesn't remove the general one.
#' @references O. Morgenstern and M. A. Woodbury, “Stability of Inverses of Input-Output Matrices,” Econometrica 18 (1950): 190.
.wb_identity <- function( A, W, I ) {
  # A: being either Lambda_s[[s]] or Phi matrix
  # I: Study (un-)specific Identity matrix
  # W: p x p matrix, or a length-p vector standing in for diag(W)
  if ( is.matrix( W ) ) {
    cp <- crossprod( A, W )    # to save computation time
    t1 <- cp %*% A
  } else {
    AW <- A * W                  # p x q, == diag(W) %*% A via row recycling
    cp <- t( AW )                 # q x p, == t(A) %*% diag(W)
    t1 <- crossprod( AW, A )      # q x q, == cp %*% A
  }
  inv_t <- solve( I + t1 )
  inv_t %*% cp
}
#' @keywords internal
#' 2nd application of the Woodbury Identity. See \code{.wb_identity()} for
#' the vector-vs-matrix \code{W} convention.
.wb_identity2 <- function( W, A, I ){
  wb <- .wb_identity( A, W, I )
  if ( is.matrix( W ) ) {
    t1 <- W %*% A
    t2 <- t1 %*% wb
    W - t2
  } else {
    t1 <- A * W                  # p x q, == diag(W) %*% A
    t2 <- t1 %*% wb               # p x p, dense result either way
    diag( W ) - t2
  }
}

#' @keywords internal
.exp_values <- function(
  Phi,              # matrix
  Lambda_s,         # list of matrices
  Psi_s,            # list of matrices
  CM_step   = 1,    # out of (1,2,3,4)
  cov_s     = NULL, # not needed for cm_step == 4
  X_s_tilde = NULL  # only needed for cm_step == 4
)
{
  k   <- dim( Phi )[2]
  I_k <- diag( 1, k )
  I_j <- lapply( Lambda_s, function( s ) diag( 1, dim( s )[2] ) )

  inv_Psi_s <- lapply( Psi_s, .inv_Psi )

  wb_f <- Map( .wb_identity2, inv_Psi_s, Lambda_s, I_j )
  wb_l <- Map( .wb_identity2, inv_Psi_s, list(Phi), list(I_k) )

  delta_Phi    <- Map( .wb_identity, list(Phi), wb_f, list(I_k) )
  delta_Lambda <- Map( .wb_identity, Lambda_s, wb_l, I_j )

  if (CM_step == 1) {
    return(list(
      exp_xl = Map( .get_exp_xl, delta_Lambda, cov_s ),
      exp_xf = Map( .get_exp_xf, delta_Phi, cov_s ),
      exp_ll = Map( .get_exp_ll, Lambda_s, delta_Lambda, cov_s, I_j ),
      exp_ff = Map( .get_exp_ff, list(Phi), delta_Phi, cov_s, list(I_k) ),
      exp_fl = Map( .get_exp_fl, Lambda_s, delta_Phi, delta_Lambda, cov_s )
    ))
  }

  if (CM_step == 2) {
    return(list(
      exp_xf = Map( .get_exp_xf, delta_Phi, cov_s ),
      exp_ff = Map( .get_exp_ff, list(Phi), delta_Phi, cov_s, list(I_k) ),
      exp_fl = Map( .get_exp_fl, Lambda_s, delta_Phi, delta_Lambda, cov_s )
    ))
  }

  if (CM_step == 3) {
    return(list(
      exp_xl = Map( .get_exp_xl, delta_Lambda, cov_s ),
      exp_ll = Map( .get_exp_ll, Lambda_s, delta_Lambda, cov_s, I_j ),
      exp_fl = Map( .get_exp_fl, Lambda_s, delta_Phi, delta_Lambda, cov_s )
    ))
  }

  if (CM_step == 4) {
    return(list(
      exp_f = Map( .get_exp_f, delta_Phi, X_s_tilde ),
      exp_l = Map( .get_exp_l, delta_Lambda, X_s_tilde )
    ))
  }
}