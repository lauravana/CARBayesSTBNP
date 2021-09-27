# y = Y.DA.mat;Xl = X.standardised; X=X.array ; w = phi.mat;
# offset = offset.mat; s = s_alloc; beta = beta.mat;
# alpha = alpha_BNP; sig2 <- nu2

sample_SFM <- function(y, Xl, X, w, offset, beta0, s, J,
                       beta, sig2, alpha, hyperparameters_P0){
  
  I <- dim(y)[1]
  N <- dim(y)[2]
  p <- dim(X)[3]
  s <- factor(s, levels = 1:J)
  nj <- unname(table(s))
  beta.offset <- (y - offset - w - beta0)

  ###########################################################
  ## Step 1: sample beta conditional on mu_beta and s
  Omega_beta <-  hyperparameters_P0$Omega_beta
  mu_all     <-  rep(hyperparameters_P0$mu_beta, J) 
  Omega_all  <-  diag(rep(diag(Omega_beta), J)) # TODO: Matrix::bdiag(rep(list(Omega_beta), J))
  s_alloc_vec <- rep(s, N)
  Xs <-  model.matrix(~ - 1 + Xl : (s_alloc_vec)) 
  S_beta <- solve(Omega_all + crossprod(Xs)/sig2)
  yX <- drop(crossprod(Xs, as.numeric(beta.offset)))
  m_beta <- (drop(mu_all %*% Omega_all) + yX/sig2) %*% S_beta
  
  beta <- MASS::mvrnorm(1, mu = m_beta, Sigma = S_beta)
  beta.mat <- matrix(beta, nrow = J, byrow = TRUE)


  
  ###########################################################
  ## Step 3: sample weights
  eta <- gtools::rdirichlet(1, alpha + nj)
  
  ###########################################################
  ## Step 3a: sample alpha
  ## assume prior: alpha ~ gamma(a, a * J)
  ## posterior is then: dgamma(alpha, a, a*J) * gammafunc(K*alpha)/gammafunc(alpha)^K*(prod(eta_k)^(alpha-1))
  ## a random walk Metropolis-Hastings step
  
  #MH step for alpha
  a <- 10
  alpha_new <- exp(log(alpha) + 0.5 * rnorm(1))
  lf <- function(x) dgamma(x, a, a * J, log = TRUE) + lgamma(x * J) - J * lgamma(x) + (x - 1) * sum(log(eta))
  log_accept <- (alpha) - (alpha_new) + lf(alpha_new) - lf(alpha)
  
  accept <- 1
  if(is.nan(log_accept)){
    accept <- 0
  }else if( log_accept < 0 ){
    accept <- exp(log_accept)
 }
  if( runif(1) < accept ){
    alpha <- alpha_new
  }
  print(alpha)
  ###########################################################
  ## Step 4: sample allocations
  probs_prop  <-  sapply(1:J, function(j) {
    dj <- dnorm(c(y), beta0 + Xl %*% beta.mat[j,] + c(w), sqrt(sig2))
    sapply(1:I, function(i) sum(dj[(1:N - 1) * I + i]))
  })
  #####
  probs  <- t(apply(probs_prop, 1, function(x) (eta * x)/sum(eta * x)))
  s <- sapply(1:nrow(probs), function(i) sample(1:J, 1, prob = probs[i, ]))
 
  # Step 5 do the random relabeling
  rho_relabel <- sample(1:J)
  
  # mu <- mu[rho_relabel, ]
  s  <- rho_relabel[s]
  
  return(list(s = s, beta = beta.mat, alpha = alpha, rho_relabel=rho_relabel))
}
