#' @keywords internal
.build_Sig <-  function( Phi, Lambda_s, Psi_s ){
  # simple generator for Sigma matrix 
  tcrossprod( Phi ) + tcrossprod( Lambda_s ) + Psi_s
}

#' @importFrom statmod vecmat 
#' @importFrom statmod matvec
#' @keywords internal
.inv_Sig <- function( Psi_s1, Lambda_s, Phi ){
  # Since direct calculation of Sigma's inverse is currelty only needed once its okay to have it as a seperate function
  # until I figure out how to retrieve it from the Woodbury Identity
  k <- dim( Phi )[2]
  j_s <- dim( Lambda_s )[2]

  LambTOT <- cbind( Lambda_s, Phi )                    # build once
  psi_vec <- diag( Psi_s1 )                            # extract once
  PsiInv_U <- statmod::vecmat( psi_vec, LambTOT )      # compute once

  inv <- solve( diag( 1, k + j_s ) + crossprod( LambTOT, PsiInv_U ) )

  Psi_s1 - PsiInv_U %*% tcrossprod( inv, PsiInv_U )
}

#'@keywords internal
.update_Psi <- function(cov, Phi, Lambda, exp_ff, exp_ll, exp_xf, exp_xl, exp_fl) {
  A <- cov + 
    Phi %*% tcrossprod( exp_ff, Phi ) + 
    Lambda %*% tcrossprod( exp_ll, Lambda ) - 
    2 * tcrossprod( exp_xf, Phi ) -  
    2 * tcrossprod( exp_xl, Lambda ) +
    2 * Phi %*% tcrossprod( exp_fl, Lambda )
  
  diag(diag(A))
}

#'@keywords internal
.update_Lambda <- function( Phi, exp_xl, exp_fl, exp_ll, j ) {
  p  <- dim(Phi)[1]
  A_ <- Phi %*% exp_fl
  B_ <- exp_xl - A_
  C_ <- solve( exp_ll )
  D_ <- B_ %*% C_

  matrix( D_, p, j )
}

#'@keywords internal
.update_first_beta_part <- function( X_og, Phi, Lambda, exp_f, exp_l, B_s ) {
  X_og_t <-  t( X_og )
  A_ <- Phi %*% exp_f
  C_ <- Lambda %*% exp_l
  D_ <- X_og_t - A_ - C_

  D_ %*% B_s
}

#' @keywords internal
.change_X <- function( X_og, B_s, beta ) {
  X_og - tcrossprod( B_s, beta )
}

#' Estimates the parameters of a MSFR model
#'
#' Maximum likelihood estimation of the MSFR model parameters via the ECM
#' algorithm.
#'
#' There are two different constraints for achieving model identification,
#' as detailed in the reference,
#' though the function can also be run without such constraints (not recommended).
#' No checking is done on the starting value for the various model matrices,
#' since a suitable value for them  is produced by the function \code{start_msfa}.
#' @param X_s List of lenght \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all the studies.
#' @param B_s List of length \eqn{S}{S}, corresponding to the number of different studies considered. 
#' Each element of the list contains a data matrix with known covariates or batch size idicators etc. 
#' @param start A list containing the slots \code{Phi}, \code{Lambda_s} and \code{Psi_s}, containing the starting
#' values for the matrix  \code{Phi} of common factor loadings, of size \eqn{P \times K}{P x K}, for
#' the matrices \code{Lambda_s} of study-specific factor loadings, a list of size \eqn{S}{S}  where each element
#' contains a matrix with \eqn{P \times J_s}{P x J_s}, and finally for the study-specific matrices of uniquenesses,
#' a list of size \eqn{S}{S}, where each element contains a vector of length \eqn{P}{P}.
#' Note that a suitable list of this kind is produced by \code{start_msfa}.
#' @param nIt Maximum number of iterations for the ECM algorithm. Default is 50000.
#' @param tol Tolerance for declaring convergence of the ECM algorithm. Default is 10^-7.
#' @param constraint  Constraint for ensuring identifiability. The default is "block_lower2", which
#' corresponds to the main proposal of De Vito et al. (2018). An alternative identification
#' strategy is triggered by  "block_lower1"; this is more restrictive but may work also with smaller
#' number of variables. Again, the latter strategy is mentioned in De Vito et al. (2018).
#' @param trace If \code{TRUE} then trace information is being printed every 1000 iterations of the ECM algorithm.
#' @return A list  containing the following components:
#' \item{\code{Phi},\code{Lambda_s}, \code{psi_s}}{the estimated model matrices.}
#' \item{loglik}{the value of the log likelihood function at the final estimates.}
#' \item{\code{AIC, BIC}}{model selection criteria at the estimate.}
#' \item{\code{npar}}{number of model parameters.}
#' \item{iter}{the number of ECM iterations performed.}
#' \item{constraint}{the identification constraint enforced.}
#' @export
#' @import robust
#' @importFrom stats cor cov factanal prcomp
#' @references De Vito, R., Bellio, R., Trippa, L. and Parmigiani, G. (2018). (2019). Multi-study Factor Analysis. Biometrics,  75, 337-346.
#' @references Pison, G., Rousseeuw, P.J., Filzmoser, P. and Croux, C. (2003). Robust factor analysis. Journal
#' of Multivariate Analysis, 84, 145-172.
ecm_msfr <- function(X_s, B_s, start, nIt = 50000, tol = 10^-7, 
                     constraint = "block_lower2", trace = TRUE)
{
  S <- length(X_s)
  ####### extract elements from input list 'start'
  Phi <- start$Phi
  Lambda_s <- start$Lambda_s
  psi_s <- start$psi_s
  beta <- start$beta

  
  # get some basic variables needed for computation
  j_s <- n_s <- numeric(S)
  p <- dim(Phi)[1]
  k <- dim(Phi)[2]
  p_b <- dim(beta)[2]
  B <- do.call( rbind, B_s )

  second_part <- solve(crossprod(B))
  theta <- .param2vect(start, constraint)

  #######defining objects
  Psi_s1 <- Psi_s <- cov_s <- list()
  L_s <- list()
  X_s_original <- X_s
  #Changing Xs for Xtilde
  X_s <- list()
  for(s in 1:S) X_s[[s]] <- X_s_original[[s]] - B_s[[s]]%*%t(beta)

  ######1st round of cycle
  for(s in 1:S){
  	n_s[s] <-  dim(X_s[[s]])[[1]]
  	j_s[s] <-  dim(Lambda_s[[s]])[[2]]
  	Psi_s[[s]] <- diag(psi_s[[s]])
    Psi_s1[[s]] <- diag(1/psi_s[[s]])
    cov_s[[s]] <- cov(X_s[[s]])
  }
  ######E-step
  # Sig_s <- Map( .build_Sig, list( Phi ), Lambda_s, Psi_s )
  # Sig_s1 <- Map( .inv_Sig, Psi_s1, Lambda_s, list(Phi) )
  # ds_s <- lapply( Sig_s, det )
  l_stop0 <- 0
  lm1 <- 0
  # l0 <- .loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
  l0 <- .loglik_int(theta, n_s, cov_s, k , j_s, constraint)

  # l.df <- data.frame(lm1 = lm1, l0 = l0, l1 = 0, a = 0, l_stop = 0) # for testing l updates
  # Algorithm loop ----
  for (i in (1:nIt)) {
    if (i%%1000 == 0){print(i)}
    ## CM1 ----
    ### expected values ----
    out <- .exp_values( Phi, Lambda_s, Psi_s, CM_step = 1, cov_s = cov_s )
    exp_xl <- out$exp_xl; exp_xf <- out$exp_xf; exp_ll <- out$exp_ll; exp_ff <- out$exp_ff; exp_fl <- out$exp_fl
    
    ### update  of Psi_s ----
    Psi_new <- Map( .update_Psi, cov_s, list(Phi), Lambda_s, exp_ff, exp_ll, exp_xf, exp_xl, exp_fl )
    Psi_new1 <- Map( .inv_Psi, Psi_new )

    ## CM2 ----
    ### expected values ----
    out <- .exp_values( Phi, Lambda_s, Psi_new, CM_step = 2, cov_s = cov_s )
    exp_xf <- out$exp_xf; exp_ff <- out$exp_ff; exp_fl <- out$exp_fl

    ### update of Phi ----
    nPsi1 <- Map( "*", n_s, Psi_new1 )
    C_s <- Map( 
      function( nPsi1, Lambda, exp_xf, exp_fl ) nPsi1 %*% ( exp_xf - tcrossprod( Lambda, exp_fl ) ),
      nPsi1, Lambda_s, exp_xf, exp_fl
    )
    kron_s <- Map(
      function( nPsi1, exp_ff ) kronecker( t(exp_ff), nPsi1 ),
      nPsi1, exp_ff
    )
    C <- Reduce( '+', C_s )
    kron <- Reduce( '+', kron_s )
    Phi_vec <- solve( kron ) %*% matrix( as.vector( C ) )
    Phi_new <- matrix( Phi_vec, p, k )

    ## CM3 ----
    ### expected values ----
    out <- .exp_values( Phi_new, Lambda_s, Psi_new, CM_step = 3, cov_s = cov_s )
    exp_xl <- out$exp_xl
    exp_ll <- out$exp_ll 
    exp_fl <- out$exp_fl

    ### update of Lambda ----
    Lambda_new <- Map(
      .update_Lambda,
      list( Phi_new ),
      exp_xl, 
      exp_fl, 
      exp_ll, 
      j_s
    )

    ## CM4 ----
    ### expected values ----
    out <- .exp_values( Phi_new, Lambda_new, Psi_new, CM_step = 4, X_s_tilde = X_s )
    exp_f <- out$exp_f; exp_l <- out$exp_l
   
    ### update of beta ----
    first_part_s <- Map(
      .update_first_beta_part,
      X_s_original,
      list(Phi_new),
      Lambda_new,
      exp_f,
      exp_l,
      B_s
    )
   
    first_part <- Reduce('+', first_part_s)
    beta_new <- first_part %*% second_part
   
    ### Changing Xs ----
    X_s <- Map(
      .change_X,
      X_s_original,
      B_s,
      list(beta_new)
    )
	  
    ## constraint ----
    lambda_vals <- c()
    psi_vals <- psi_new <- c()
    Phi_new[upper.tri(Phi_new)] <- 0
    phi_val <- as.vector(Phi_new[lower.tri(Phi_new, diag = TRUE)])
    for (s in 1:S){
      Lambda_new[[s]][upper.tri(Lambda_new[[s]])] <- 0
      lambda_vals <- c(lambda_vals, as.vector(Lambda_new[[s]][lower.tri(Lambda_new[[s]], diag = TRUE)]))
      psi_new[[s]] <- diag(Psi_new[[s]])
      psi_vals <- c(psi_vals, psi_new[[s]])
    }
    L_sTOT <- Reduce('cbind', Lambda_new)
    Omega <- cbind(Phi_new, L_sTOT)
    rank_tot <-  qr(Omega)$rank
    theta_new <- theta_test <- c(phi_val, lambda_vals, psi_vals)
    param.struct <- list(Phi = Phi_new, Lambda_s = Lambda_new, psi_s=psi_new)
    Delta <- theta_new - theta
    sh <- 0   ###no more than 20 step-halving rounds

    while( (rank_tot < k + sum(j_s)) & (sh<20)) {
      Delta <- Delta / 2
      sh <- sh + 1
      theta_test <- theta + Delta
      param <- .vect2param(theta_test, param.struct, constraint, p, k, j_s)
      Lambda_new <- c()
      psi_new <- param$psi_s
      for(s in 1:S) {
        Lambda_new[[s]] <- param$Lambda_s[[s]]
        Psi_new[[s]] <- diag(psi_new[[s]])
        Psi_new1[[s]] <- diag(1 / psi_new[[s]])
      }
      L_sTOT <- Reduce('cbind', Lambda_new)
      Phi_new <- param$Phi
      Omega <- cbind(Phi_new, L_sTOT)
      rank_tot <-  qr(Omega)$rank
    }

    if(sh==20) stop("The full rank condition does not hold\n")


    ## stopping rule ----
    # Sig_s <- Map( .build_Sig, list( Phi ), Lambda_s, Psi_s )
    # Sig_s1 <- Map( .inv_Sig, Psi_new1, Lambda_new, list(Phi_new) ) 
    # ds_s <- lapply( Sig_s, det )

    # l1 <- .loglik_ecm(Sig_s1, ds_s, n_s, cov_s)
    l1 <- .loglik_int(theta_new, n_s, cov_s, k , j_s, constraint)
    a <- (l1 - l0)/ (l0-lm1)
    l_stop <- lm1 + (1/ (1-a)) * (l0-lm1)
    
    # save stopping values for debugging
    # df <-  data.frame(lm1 = lm1, l0 = l0, l1 = l1, a = a, l_stop = l_stop)
    # l.df <- dplyr::bind_rows(l.df, df)

    # l0 <- .loglik_ecm(Sig_s1,  ds_s, n_s, cov_s) # seemed odd 
    if((trace) & (i %% 1000 == 0))  cat("i=", i, "Criterion for convergence ", abs(l_stop-l_stop0),  "\n")
    if((abs(l_stop-l_stop0)<tol) & i > 1 & l_stop != Inf) break

    # assign vars for new cycle
    Psi_s <- Psi_new
    psi_s <- psi_new
    Phi <- Phi_new
    Lambda_s <- Lambda_new
    beta <- beta_new
    if (constraint == "block_lower2") theta <- theta_new
    Psi_s1 <- Psi_new1
    lm1 <- l0
    l0 <- l1
    l_stop0 <- l_stop
  }
  
  # AIC and BIC computation ----

  if (constraint == "block_lower1") npar <- p * S + k * (p - ( k - 1) / 2) +  sum(j_s * (p - k - (j_s - 1) / 2))
  if (constraint == "block_lower2")  npar <- p * S + k * (p - ( k - 1) / 2) +  sum(j_s * (p  - (j_s - 1) / 2))
  n_tot <- sum(n_s)
  AIC <- -2 * l1 + npar * 2
  BIC <- -2 * l1 + npar * log(n_tot)

  # return output ----
  res <- list(Phi = Phi, Lambda_s = Lambda_s, beta = beta, psi_s = psi_s, loglik = l1,
              AIC = AIC, BIC = BIC, npar=npar,
              iter = i,  cov_s = cov_s,  n_s = n_s, constraint=constraint)
  # write.csv2(l.df, "likelihoods_no253.csv")
  return(res)
}
