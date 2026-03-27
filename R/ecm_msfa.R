#' @keywords internal
#' simple generator for Sigma matrix 
.build_Sig <-  function( Phi, Lambda_s, Psi_s ){
  tcrossprod( Phi ) + tcrossprod( Lambda_s ) + Psi_s
}

#' @importFrom statmod vecmat 
#' @importFrom statmod matvec
#' @keyword internal
#' Since direct calculation of Sigma's inverse is currelty only needed once its okay to have it as a seperate function
#' until I figure out how to retrieve it from the Woodbury Identity
.inv_Sig <- function( Psi_s1, Lambda_s, Phi ){
  k <- dim( Phi )[2]
  S <- length( Lambda_s )
  j_s <-  sapply( Lambda_s, function( l ) dim( l )[2] )
  I_tot <- lapply( k + j_s, diag, x = 1 )
  LambTOT <-  lapply( Lambda_s, cbind, Phi )
  list_of_args <- .make_args( list( Psi_s1, I_tot, LambTOT ) )

  lapply(
    list_of_args,
    do.call, 
    what = function( ps1, i_tot, lamtot ) {
      ps1 - (
        statmod::vecmat( diag( ps1 ), lamtot ) %*%
          solve(
            i_tot + (
              t( lamtot ) %*% 
                statmod::vecmat( diag( ps1 ), lamtot )
              )
          ) %*% 
        statmod::matvec( t( lamtot ), diag( ps1 ) )
      )             
    }
  ) 
}


#' Estimates the parameters of a MSFA model
#'
#' Maximum likelihood estimation of the MSFA model parameters via the ECM
#' algorithm.
#'
#' There are two different constraints for achieving model identification,
#' as detailed in the reference,
#' though the function can also be run without such constraints (not recommended).
#' No checking is done on the starting value for the various model matrices,
#' since a suitable value for them  is produced by the function \code{start_msfa}.
#' @param X_s List of lenght \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all the studies.
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
#' @param robust If \code{TRUE}, robust covariance matrix is used in place of the sample covariance. Default
#' is \code{FALSE}.
#' @param corr If \code{TRUE}, the analysis will employ the correlation matrix instead of the covariance matrix.
#' @param mcd If \code{TRUE}, the robust estimator used for the covariance is the same proposed in Pison et al. (2003),
#' otherwise the default value of the function \code{CovRob} of the \code{robust} library is employed. Default is
#' \code{FALSE}.
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
ecm_msfa <- function(X_s, B_s, start, nIt = 50000, tol = 10^-7, 
                     constraint = "block_lower2", trace = TRUE)
{
  S <- length(X_s)
  ####### extract elements from input list 'start'
  Phi <- replicate( S, start$Phi, simplify = F )
  Lambda_s <- start$Lambda_s
  psi_s <- start$psi_s
  beta <- replicate( S, start$beta, simplify = F )

  
  # get some basic variables needed for computation
  j_s <- n_s <- numeric(S)
  p <- dim(Phi[[1]])[1]
  k <- dim(Phi[[1]])[2]
  p_b <- dim(beta[[1]])[2]
  B <- do.call( rbind, B_s )

  second_part <- solve(t(B) %*% B)
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
    Psi_s1[[s]] <-  diag(1/psi_s[[s]])
    cov_s[[s]] <- cov(X_s[[s]])
  }
  ######E-step
  Sig_s <- Map( .build_Sig, list( Phi ), Lambda_s, Psi_s )
  Sig_s1 <- .inv_Sig( Psi_s1, Lambda_s, Phi )
  ds_s <- lapply( Sig_s, det )
  l_stop0 <- 0
  lm1 <- 0
  l0 <- .loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
  for (i in (1:nIt)) {
    if (i%%100 == 0){print(i)}
    ###########CM1 ---------------------------------------------------------------------------------------

    ###### expected values
    out <- .exp_values( Phi, Lambda_s, Psi_s, CM_step = 1, cov_s = cov_s )
    exp_xl <- out$exp_xl; exp_xf <- out$exp_xf; exp_ll <- out$exp_ll; exp_ff <- out$exp_ff; exp_fl <- out$exp_fl
    ###### update  of Psi_s
    Psi_new <- list()
    Psi_new1 <- list()
    psi_new <- list()


    for(s in 1:S){
   	  psi_new[[s]]  <- diag(cov_s[[s]] + Phi %*% exp_ff[[s]] %*% t(Phi) + Lambda_s[[s]] %*%
   	                   exp_ll[[s]] %*% t(Lambda_s[[s]]) - 2*exp_xf[[s]] %*% t(Phi) -  2*exp_xl[[s]] %*% t(Lambda_s[[s]]) +
   	                   2 * Phi %*% exp_fl[[s]] %*% t(Lambda_s[[s]]))
   	  Psi_new[[s]] <- diag(psi_new[[s]])
   	  ##########inverse
   	  Psi_new1[[s]] <- diag(1/diag(Psi_new[[s]]))
   	}

    ###########CM2 ---------------------------------------------------------------------------------------

    ######expected values
    out <- .exp_values(Phi, Lambda_s, Psi_new, CM_step = 2, cov_s = cov_s )
    exp_xl <- out$exp_xl; exp_xf <- out$exp_xf; exp_ll <- out$exp_ll;
    exp_ff <- out$exp_ff; exp_fl <- out$exp_fl

    ######update of Phi
    C_s <- list()
    kron_s <- list()
    for(s in 1:S){
      C_s[[s]] <- n_s[s] * Psi_new1[[s]] %*% exp_xf[[s]] - n_s[s] * Psi_new1[[s]] %*% Lambda_s[[s]] %*% t(exp_fl[[s]])
      kron_s[[s]] <- kronecker(t(exp_ff[[s]]), n_s[s] * Psi_new1[[s]])
    }
    C <- Reduce('+', C_s)
    kron <- Reduce('+', kron_s)
    Phi_vec <- solve(kron) %*% matrix(as.vector(C))
    Phi_new <- matrix(Phi_vec, p, k)


    ########CM3 ---------------------------------------------------------------------------------------

    ######expected values
    out <- .exp_values( Phi_new, Lambda_s, Psi_new, CM_step = 3, cov_s = cov_s )
    exp_xl <- out$exp_xl; exp_xf <- out$exp_xf; exp_ll <- out$exp_ll;
    exp_ff <-  out$exp_ff; exp_fl <- out$exp_fl

    ######update of Lambda
    Lambda_new <- list()
    for(s in 1:S){
      Lambda_new[[s]] <- matrix(((exp_xl[[s]] - Phi_new %*% exp_fl[[s]]) %*% solve(exp_ll[[s]])), p, j_s[s])
    }

    ########CM4: new part for beta---------------------------------------------------------------------------------------
   
    ######expected values
    out <-  .exp_values( Phi_new, Lambda_new, Psi_new, CM_step = 4, X_s_tilde = X_s )
    exp_xl <- out$exp_xl; exp_xf <- out$exp_xf; exp_ll <- out$exp_ll;
    exp_ff <-  out$exp_ff; exp_fl <- out$exp_fl; exp_f <- out$exp_f; exp_l <- out$exp_l
   
    ######update of beta
    first_part_s <- list()
    for (s in 1:S){
     first_part_s[[s]] <- (t(X_s_original[[s]]) - Phi_new %*%  exp_f[[s]] - Lambda_new[[s]] %*% exp_l[[s]]) %*% B_s[[s]]		
    }
   
    first_part <- Reduce('+', first_part_s)
    beta_new = first_part %*% second_part
   
    #Changing Xs
    X_s <- list()
    for(s in 1:S) X_s[[s]] <- X_s_original[[s]] - B_s[[s]]%*%t(beta_new)
	  
    ###### constraint ---------------------------------------------------------------------------------------------------
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
    theta_new <- c(phi_val, lambda_vals, psi_vals)
    param.struct <- list(Phi = Phi_new, Lambda_s = Lambda_new, psi_s=psi_new)
    Delta <- theta_new - theta
    sh <- 0   ###no more than 20 step-halving rounds

    while( (rank_tot < k + sum(j_s)) & (sh<20)) {
      Delta <- Delta / 2
      sh <- sh + 1
      theta_new <- theta + Delta
      param <- .vect2param(theta_new, param.struct, constraint, p, k, j_s)
      Lambda_new <- c()
      psi_new <- param$psi_new
      for(s in 1:S) {
        Lambda_new[[s]] <- param$Lambda_s[[s]]
        Psi_new[[s]] <- diag(psi_new[[s]])
        Psi1_new[[s]] <- diag(1 / psi_new[[s]])
      }
      L_sTOT <- Reduce('cbind', Lambda_new)
      Phi_new <- param$Phi
      Omega <- cbind(Phi_new, L_sTOT)
      rank_tot <-  qr(Omega)$rank
    }

    if(sh==20) stop("The full rank condition does not hold\n")


    ###### stopping rule -------------------------------------------------------------------------------------------
    Sig_s <- Map( .build_Sig, list( Phi ), Lambda_s, Psi_s )
    Sig_s1 <- .inv_Sig( Psi_new1,Lambda_new, Phi_new ) 
    ds_s <- lapply( Sig_s, det )

    l1 <- .loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
    a <- (l1 - l0)/ (l0-lm1)
    l_stop <- lm1 + (1/ (1-a)) * (l0-lm1)
    
    l0 <- .loglik_ecm(Sig_s1,  ds_s, n_s, cov_s)
    if((trace) & (i %% 1000 == 0))  cat("i=", i, "Criterion for convergence ", abs(l_stop-l_stop0),  "\n")
    if((abs(l_stop-l_stop0)<tol) & i > 1 & l_stop != Inf) break
    
    # assign vars for new cycle
    Psi_s <- Psi_new
    Phi <- Phi_new
    Lambda_s <- Lambda_new
    beta <- beta_new
    if (constraint == "block_lower2") theta <- theta_new
    Psi_s1 <- Psi_new1
    lm1 <- l0
    l0 <- l1
    l_stop0 <- l_stop
  }
  #####AIC and BIC computation

  if (constraint == "block_lower1") npar <- p * S + k * (p - ( k - 1) / 2) +  sum(j_s * (p - k - (j_s - 1) / 2))
  if (constraint == "block_lower2")  npar <- p * S + k * (p - ( k - 1) / 2) +  sum(j_s * (p  - (j_s - 1) / 2))
  n_tot <- sum(n_s)
  AIC <- -2 * l1 + npar * 2
  BIC <- -2 * l1 + npar * log(n_tot)

  ############return output
  res <- list(Phi = Phi, Lambda_s = Lambda_s, beta = beta, psi_s = psi_s, loglik = l1,
              AIC = AIC, BIC = BIC, npar=npar,
              iter = i,  cov_s = cov_s,  n_s = n_s, constraint=constraint)
  return(res)
}
