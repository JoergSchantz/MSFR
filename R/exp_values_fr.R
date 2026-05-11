
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

    Sig_s[[s]] <- tcrossprod( Phi ) + Psi_s[[s]]

    if ( getdet ) ds_s[[s]] <- det( Sig_s[[s]] )

    Sig_s1[[s]] <- Psi_s1[[s]] - (statmod::vecmat(diag(Psi_s1[[s]]), Phi) %*%
                                    solve(I_k + (t(Phi) %*% statmod::vecmat(diag(Psi_s1[[s]]),
                                                                                             Phi ))) %*% statmod::matvec(t( Phi ), diag(Psi_s1[[s]])))

    delta_Phi[[s]] <- crossprod( Phi, Sig_s1[[s]] )

    Delta_Phi[[s]] <- I_k - ( crossprod( Phi, Sig_s1[[s]] ) %*% Phi )

    Txsfcs[[s]] <- tcrossprod( cov_s[[s]], delta_Phi[[s]] )

    Tfcsfcs[[s]] <- delta_Phi[[s]] %*% tcrossprod( cov_s[[s]], delta_Phi[[s]] ) + Delta_Phi[[s]]

    ##new part for beta

    Woodbury_f[[s]] <- .wb_identity( Phi, .inv_Psi( Psi_s[[s]] ), I_k )

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