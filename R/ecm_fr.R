#' Estimates the parameters of study-specific FR models
#'
#' Maximum likelihood estimation of study-specific FR models parameters via the ECM
#' algorithm, adopting the upper-triangular zero constraint to achieve identification
#' for each loading matrix. Note: the function can also estimate a FR model for a single
#' study, by specifiyng \code{X_s = list(data)}, where \code{data} is the data matrix.
#' @param X_s List of lenght \eqn{S}{S}, corresponding to number of different studies considered.
#' Each element of the list contains a data matrix, with the same number of columns \eqn{P}{P} for all the studies.
#' @param B_s List of length \eqn{S}{S}, corresponding to the number of different studies considered. 
#' Each element of the list contains a data matrix with known covariates or batch size idicators etc. 
#' @param start initialised parameters i.e. output from \code{start_fr()}
#' @param nIt Maximum number of iterations for the ECM algorithm. Default is 50000.
#' @param tol Tolerance for declaring convergence of the ECM algorithm. Default is 10^-7.
#' @param block_lower Should the upper-triangular zero constraint be enforced? Default is \code{TRUE}
#' (strongly suggested).
#' @param robust If \code{TRUE}, robust covariance matrix is used in place of the sample covariance. Default
#' is \code{FALSE}.
#' @param corr If \code{TRUE}, the analysis will employ the correlation matrix instead of the covariance matrix.
#' @param mcd If \code{TRUE}, the robust estimator used for the covariance is the same proposed in Pison et al. (2003),
#' otherwise the default value of the function \code{CovRob} of the \code{robust} library is employed. Default is
#' \code{FALSE}.
#' @param trace If \code{TRUE} then trace information is being printed every \code{traceIT} iterations of the ECM algorithm.
#' @import robust psych
#' @references De Vito, R., Bellio, R., Trippa, L. and Parmigiani, G. (2019). Multi-study Factor Analysis. Biometrics,  75, 337-346.
#' @references Pison, G., Rousseeuw, P.J., Filzmoser, P. and Croux, C. (2003). Robust factor analysis. Journal
#' Multivariate Analysis, 84, 145-172.
#' @return A list  containing the following components:
#' \item{\code{Phi_s}, \code{beta} and \code{Psi_s}}{the estimated model matrices.}
#' \item{loglik}{the value of the log likelihood function at the final estimates.}
#' \item{\code{AIC, BIC}}{model selection criteria at the estimate.}
#' \item{\code{npar}}{number of model parameters.}
#' \item{iter}{the number of ECM iterations performed.}
#' @examples 
#' data(Scenario1_MSFR)
#' EM <- ecm_fr(X_s, B_S)
#' # estimated study-specific model matrices
#' Phi_1 <- EM$Phi_s[[1]]
#' # ...
#' # regression matrix
#' beta <- EM$beta
#' # study-specific error matrices
#' Psi_1 <- EM$Psi_s[[1]]
#' # ...
#' # visualise these matrices (here only for one)
#' heat_plot(Phi_1)
#' @export
ecm_fr <- function( 
  X_s,
  B_s,
  start,
  nIt = 50000,
  tol = 10 ^ -7,
  block_lower = TRUE,
  robust = FALSE,
  corr = TRUE,
  mcd = FALSE,
  trace = TRUE
)
{
  Psi_s <- psi_s <- list()
  #######
  p     <- ncol( X_s[[1]] )
  k     <- dim( start$Phi )[[2]]
  S     <- length( X_s )
  n_s   <- numeric( S )
  #######defining objects
  Psi_s1 <- Psi_s <- cov_s <- list()

  tot_s <- k

  Phi   <- start$Phi
  beta  <- start$beta
  psi_s <- start$psi_s

  B           <- Reduce( 'rbind', B_s )
  second_part <- solve( crossprod( B ) )

  X_s_original <- X_s
  #Changing Xs for Xtilde
  X_s          <- list()
  for( s in 1:S ) X_s[[s]] <- X_s_original[[s]] - B_s[[s]] %*% t( beta )

  ######1st round of cycle
  for( s in 1:S ) {
    n_s[s]       <- dim( X_s[[s]] )[[1]]
    Psi_s[[s]]   <- diag( psi_s[[s]] )
    Psi_s1[[s]]  <- .inv_Psi( Psi_s[[s]] )
    cov_s[[s]]   <- cov( X_s[[s]] )
  }
  ######E-step
  out    <- .exp_values_fr(
              Phi,
              Psi_s,
              cov_s,
              X_s,
              getdet = TRUE
            )
  Sig_s1  <- out$Sig_s1
  ds_s    <- out$ds_s
  l_stop0 <- 0
  lm1     <- 0
  l0      <- .loglik_ecm( Sig_s1, ds_s, n_s, cov_s )
  for( i in ( 1:nIt ) )
  {
    ###########CM1 ---------------------------------------------------------------------------------------

    ######expected values
    out     <- .exp_values_fr( Phi,
                                Psi_s,
                                cov_s,
                                X_s )
    Txsfcs  <- out$Txsfcs
    Tfcsfcs <- out$Tfcsfcs
    ######update  of Phi_s
    Psi_new  <- list()
    Psi_new1 <- list()
    psi_new  <- list()

    for( s in 1:S ) {
      term2         <- Phi %*% Tfcsfcs[[s]] %*% t( Phi )
      term3         <- 2 * Txsfcs[[s]] %*% t( Phi )
      psi_new[[s]]  <- cov_s[[s]] + term2 - term3
      Psi_new[[s]]  <- diag( diag( psi_new[[s]] ) )
      ##########inverse
      Psi_new1[[s]] <- .inv_Psi( Psi_new[[s]] )
    }

    ###########CM2 ---------------------------------------------------------------------------------------

    ######expected values
    out     <- .exp_values_fr( Phi,
                                Psi_new,
                                cov_s,
                                X_s )
    Txsfcs  <- out$Txsfcs
    Tfcsfcs <- out$Tfcsfcs

    ######update of Phi
    C_s    <- list()
    kron_s <- list()
    for( s in 1:S ) {
      C_s[[s]]    <- n_s[s] * Psi_new1[[s]] %*% Txsfcs[[s]]
      kron_s[[s]] <- kronecker( t( Tfcsfcs[[s]] ), n_s[s] * Psi_new1[[s]] )
    }
    C       <- Reduce( '+', C_s )
    kron    <- Reduce( '+', kron_s )
    Phi_vec <- solve( kron ) %*% matrix( as.vector( C ) )
    Phi_new <- matrix( Phi_vec, p, k )

    ########CM3: new part for beta---------------------------------------------------------------------------------------

    ######expected values
    out        <- .exp_values_fr( Phi_new,
                                   Psi_new,
                                   cov_s,
                                   X_s )
    Txsfcs     <- out$Txsfcs
    Tfcsfcs    <- out$Tfcsfcs
    E_fis_x_is <- out$E_fis_x_is

    ######update of beta
    first_part_s <- list()
    for( s in 1:S ) {
      cross_term        <- Phi_new %*% E_fis_x_is[[s]]
      first_part_s[[s]] <- ( t( X_s_original[[s]] ) - cross_term ) %*% B_s[[s]]
    }

    first_part <- Reduce( '+', first_part_s )
    beta_new   <- first_part %*% second_part

    #Changing Xs
    X_s <- list()
    for( s in 1:S ) X_s[[s]] <- X_s_original[[s]] - B_s[[s]] %*% t( beta_new )

    ###########stopping rule
    out <- .exp_values_fr( Phi_new,
                            Psi_new,
                            cov_s,
                            X_s,
                            getdet = TRUE )

    Sig_s1  <- out$Sig_s1
    ds_s    <- out$ds_s
    l1      <- .loglik_ecm( Sig_s1, ds_s, n_s, cov_s )
    a       <- ( l1 - l0 ) / ( l0 - lm1 )
    l_stop  <- lm1 + ( 1 / ( 1 - a ) ) * ( l0 - lm1 )
    if( ( trace ) & ( i %% 100 == 0 ) ) {
      cat( "i=",
           i,
           "Criterion for convergence ",
           abs( l_stop - l_stop0 ),
           "\n" )
    }
    if( ( abs( l_stop - l_stop0 ) < tol ) & i > 1 & l_stop != Inf ) break
    Psi_s   <- Psi_new
    Phi     <- Phi_new
    Phi_s   <- Phi_new
    Beta_s  <- beta_new
    Psi_s1  <- Psi_new1
    lm1     <- l0
    l0      <- l1
    l_stop0 <- l_stop
  }
  ############return output
  for( s in 1:S ) psi_s[[s]] <- diag( Psi_s[[s]] )
  npar  <- p * S + sum( tot_s * ( p - ( tot_s - 1 ) / 2 ) )
  n_tot <- sum( n_s )
  AIC   <- -2 * l1 + npar * 2
  BIC   <- -2 * l1 + npar * log( n_tot )
  res   <- list( Phi = Phi_s,
                 Psi_s = psi_s,
                 beta = Beta_s,
                 loglik = l1,
                 AIC = AIC,
                 BIC = BIC,
                 npar = npar,
                 iter = i )
  return( res )
}
