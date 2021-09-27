#' Build precision matrix Q(rho, P) for the GMRF
#' 
#' @param n1 number of rows in grid
#' @param n2 number of columns in grid
#' @param rho spatial autocorrelation parameter between 0 and 1
#' @export
build_Q_sparse <- function(n1, n2, rho){
  diags <- list(rep(c(rep(1, (n1 - 1)), 0), n2)[-(n1 * n2)], rep(1, n1 * (n2 - 1)), 
                rep(c(rep(1, (n1 - 1)), 0), (n2-1))[- n1 * (n2 - 1)], 
                c(rep(c(0, rep(1, (n1 - 1))), (n2 - 1)),0))
  P     <- bandSparse(n1 * n2, k = -c(1, n1, (n1 + 1), (n1 - 1)), 
                      diag = diags, symm = TRUE)
  diagI <- bandSparse(n1 * n2, k = 0, diag = list(rep(1, n1 * n2)), symm=TRUE)
  sumP  <- bandSparse(n1*n2, k = 0, diag = list(as.vector(P %*% rep(1, n1 * n2))), symm=TRUE)
  
  Q     <- rho * (sumP - P) + (1 - rho) * diagI
  return(Q = Q)
}
