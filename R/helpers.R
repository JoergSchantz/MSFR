# This is Joerg's playground for the MSFR package 
## Probably will result in a helper file, which contains all 'internal' functionalities
## lets goooo

#' @keywords internal
#' since we are storing our e.g. study specific matrices, in lists its good to have
#' function that converts all of these lists into an new list of arguemnts for any 
#' usecase. This enables us to make use of the vectorization of R.
.make_args <-  function( LIST, ... ) {
  # check if all sublists have same length
  if ( ! length( unique( lengths( LIST ) ) ) == 1 ) {
    stop("Unequal number of elements!") # Improve message in the future
  }
  n <-  unique( lengths( LIST ) )
  lapply( 1:n, function( s ) { lapply( LIST, `[[`, s ) } )
}

#' @keywords internal
#' Computes trace of a matrix
.tr <- function( A ) sum( diag( A ) )

#' @keywords internal
#' Inverts often used Psi matrix more efficiently than solve()
.inv_Psi <- function( Psi ) diag( 1 / diag( Psi ) )

#' @keywords internal
#' Same as .inv_Psi(), but returns the length-p diagonal as a plain vector
#' instead of a dense p x p matrix. Psi_s^-1 is always diagonal, so callers
#' that only ever feed it into .wb_identity()/.wb_identity2() (as opposed to
#' needing a real p x p object, e.g. for kronecker() in ecm_fr.R) should use
#' this instead: it lets those functions take their vector-aware fast path
#' and skip an O(p^2) allocation entirely. See .wb_identity()'s docs.
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
.loglik_int <- function(theta, n_s, cov_s, k, j_s, constraint) {
 
  S  <- length(n_s)
  p  <- ncol(cov_s[[1]])
  nP <- k * p - k * (k - 1) / 2
 
  if (constraint == "block_lower1") nL <- j_s * (p - k) - j_s * (j_s - 1) / 2
  if (constraint == "block_lower2") nL <- j_s * p - j_s * (j_s - 1) / 2
 
  phi_vals <- theta[1:nP]
  Phi <- matrix(0, p, k)
  Phi[lower.tri(Phi, diag = TRUE)] <- phi_vals
 
  ## index bookkeeping done once, not inside the loop (was O(S^2), now O(S))
  nL_cum   <- cumsum(c(0, nL))     # nL_cum[s] = offset before group s
  totalL   <- nL_cum[S + 1]
  psi_base <- nP + totalL
 
  group_contrib <- function(s) {
 
    ind   <- (nP + nL_cum[s] + 1):(nP + nL_cum[s] + nL[s])
    ind_s <- (psi_base + p * (s - 1) + 1):(psi_base + p * s)
    Omega_s <- matrix(0, p, k + j_s[s])
 
    if (constraint == "block_lower1") {
      omega_vals_s <- c(phi_vals, theta[ind])
      Omega_s[lower.tri(Omega_s, diag = TRUE)] <- omega_vals_s
    }
    if (constraint == "block_lower2") { 
      Lambda_s <- matrix(0, p, j_s[s])
      Lambda_s[lower.tri(Lambda_s, diag = TRUE)] <- theta[ind]
      Omega_s <- cbind(Phi, Lambda_s)
    }
 
    psi_vals_s <- theta[ind_s]
    inv_psi_s  <- 1 / psi_vals_s

    D1L_s <- Omega_s * sqrt(inv_psi_s)
    D2L_s <- Omega_s * inv_psi_s
 
    LDL_s <- crossprod(D1L_s)                # (k+j_s) x (k+j_s)
    A     <- diag(k + j_s[s]) + LDL_s
    cholA <- chol(A)
    A1    <- chol2inv(cholA)
 
    log_ds_s <- 2 * sum(log(diag(cholA))) + sum(log(psi_vals_s))
 
    Cov_s <- cov_s[[s]]
 
    ## Woodbury trace identity: never forms the dense p x p Sig_s1
    trace_diag <- sum(inv_psi_s * diag(Cov_s))            # O(p)
    M          <- crossprod(D2L_s, Cov_s %*% D2L_s)       # O(p^2 * (k+j_s))
    trace_corr <- sum(A1 * M)                             # both symmetric
    trace_term <- trace_diag - trace_corr
 
    -(n_s[s] / 2) * log_ds_s - (n_s[s] / 2) * trace_term
  }
 
  sum(vapply(seq_len(S), group_contrib, numeric(1)))
}