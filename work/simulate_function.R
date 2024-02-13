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
  return(list(W = P, Q = Q))
}

#' Simulate data from the spatio-temporal GMRF
#' 
#' @param n1 number of rows in grid
#' @param n2 number of columns in grid
#' @param N number of time points
#' @param J number of clusters
#' @param K number of covariates to be specified if harmonic = FALSE (including intercept)
#' @param betas a JxK matrix of cluster specific regression coefficients, 
#' @param xi cluster specific temporal autoregressive coefficients,
#' @param rho spatial auto-correlation parameter between 0 and 1
#' @param tau2 variance of the spatial component
#' @param sig2 variance of the observations
#' @param tau2.delta variance of the area effects.
#' @param harmonic boolean indicating whether the covariates are 
#'   harmonic = TRUE or normal harmonic = FALSE
#' @param freq frequencies for the harmonic covariates if harmonic = TRUE
#' @param seed the seed used in the simulation
#' @param long boolean indicating whether a long data set should be provided

simulate_spatio_temporal_model <- function(n1, n2, N, J, K = NULL,
                                           betas = NULL, xi = NULL,
                                           rho, tau2, sig2, 
                                           harmonic = TRUE, 
                                           sampleBeta0BNP = TRUE,
                                           freq = NULL, 
                                           tau2.delta = 0,
                                           seed = 123,
                                           long = TRUE) {
  library(Matrix)
  library(mvtnorm)
  library(slam)
  library(spdep)
  set.seed(seed)
  
  ##################################################
  ######## FIXING THE PARAMETERS & CHECKS ##########
  ##################################################
  I <- n1 * n2 ## number of areal units in the grid

  if (!is.null(betas) && nrow(betas) != J) stop("Number of rows in betas does not match number of clusters J.")
  if (is.null(K) & !harmonic) stop("Number of normal covariates missing.")
  
  # Build regression covariates
  if (harmonic) {
    if (is.null(freq)) freq <- 7 * 1:2
    omega <- 2 * pi * freq/N ## periodicity of time series
    H <- lapply(omega, function(x)
      cbind(cos(x * seq_len(N)),
            sin(x * seq_len(N))))
    H <- do.call("cbind", H)
  } else {
    H <- cbind(matrix(rnorm(N * K), ncol = K))  ### simulated generic covariates
  }
  if (sampleBeta0BNP) H <- cbind(1, H)
  K <- ncol(H)
  

  if (!is.null(xi)) {
    if(length(xi) == 1)  xi <- rep(xi, J)
  } else {
    xi <- round(runif(J, 0.1, 0.9), 2)
  }
  
  if (is.null(betas)) {
    betas <- round(matrix(rnorm(K * J), nrow = J), 2)
  } 
  ######################
  #### Simulation ######
  ######################
  # 1) Build the sparse Q matrix
  WQ <- build_Q_sparse(n1, n2, rho = rho)
  

  # 3) sample an allocation vector of length I, which takes values in 1:J
  s <- sample(seq_len(J), I, replace = TRUE, p = 1:J/sum(1:J))
  xi.tot <- xi[s]
  betas.tot <- betas[s, ]
  beta0 <- if (sampleBeta0BNP) 0 else 3
  
  # 4) compute the expected matrix
  Hbeta <- beta0 + tcrossprod(betas.tot[, 1:ncol(H)], H)
  
  # 5) now we need to sample the spatio temporal random effects w_it
  delta <- rnorm(I, 0, sqrt(tau2.delta))
  delta.mat <- matrix(rep(delta, N), ncol = N)
  random.effects.matrix <- matrix(NA, I, N)
  epsilon.matrix <- matrix(rnorm(I * N), nrow = I)
  cholQtau <- chol(WQ$Q/tau2)
  w.demeaned <- backsolve(cholQtau, epsilon.matrix)
  random.effects.matrix[,1] <- w.demeaned[, 1]
  for (t in 2:N) {
    random.effects.matrix[, t] <- xi.tot * random.effects.matrix[, t - 1] + w.demeaned[, t]
  }
  
  # 6) sample the variance of the observations
  y.errors <- matrix(rnorm(I * N, 0, sqrt(sig2)), nrow = I)
  
  # 7) put it all together for the observations
  y.matrix <- Hbeta + random.effects.matrix + delta.mat + y.errors
  
  ## RETURN
  ## Parameters
  params <-  list(s = s, betas = betas, xi = xi, W = as.matrix(WQ$W))
  
  ### make data in long format
  if (long) {
    if (sampleBeta0BNP) H  <-  H[,-1]
    Xmat <- do.call("rbind", rep(list(H), I))
    t_id <- rep(1:N, I)
    df <- data.frame(y = c(y.matrix), Xmat[order(t_id), ],
                     area_id = rep(1:I, N), time_id = t_id[order(t_id)])
    out <- list(data = df, params = params)
  } else {
    out <- list(data = list(y.matrix, H),  params = params)
  }
  return(out)
}
