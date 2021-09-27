# y = Y.DA.mat;Xl = X.standardised; X=X.array ; w = phi.mat;
# offset = offset.mat; s = s_alloc; beta = beta.mat;
# alpha = alpha_BNP; sig2 <- nu2
# xi <- gamma

sample_BNP <- function(y, Xl, X, w, offset, beta0, s, 
                       beta, sig2, tau2, xi, rho,
                       alpha, N_aux,
                       hyperparameters_P0){

  I <- dim(y)[1]
  N <- dim(y)[2]
  p <- dim(X)[3]
  J <- length(unique(s))
  print(J)
  nj <- unname(table(s))
  # y_w_X <- lapply(1:I, function(i) (beta.offset[i, ])%*% X[i, ,])
  if (!is.null(xi)) {
    Q <- build_Q_sparse(sqrt(I), sqrt(I), rho = rho)
    w_t_1 <- w[, - ncol(w)] # this is the w_{t-1} matrix
    w_demeaned <- w[, -1] - xi[s] * w_t_1 
    s2_c <- 1/diag(as.matrix(Q))
  }

  mu_beta    <- hyperparameters_P0$mu_beta
  Omega_beta <- hyperparameters_P0$Omega_beta
  Sigma_beta <- hyperparameters_P0$Sigma_beta
  a_xi <- hyperparameters_P0$a_xi
  b_xi <- hyperparameters_P0$b_xi
  s_xi <- hyperparameters_P0$s_xi
  
  #######################
  ## Neal Algorithm 8: ##
  #######################
  # Step 1: First sample the new allocations s in a loop
  # Renovate auxiliary variables
  beta_aux <-  mvnfast::rmvn(N_aux, mu_beta, Sigma_beta)
  #beta_aux <-  matrix(rnorm(p * N_aux, rep(mu_beta, N_aux), 
  #                          rep(diag(Sigma_beta), N_aux)),
  #                    nrow=N_aux,byrow=TRUE)
  if (!is.null(xi)) xi_aux <- 2 * rbeta(N_aux, a_xi, b_xi) - 1
  for (i in seq_len(I)) {
    #Allocation of i-th location
    aux_s <- s[i]
    #Remove element from count
    nj[aux_s] <- nj[aux_s] - 1
    #It could be the only element with j-th label
    alone_j <- nj[aux_s] == 0
    #Part of the re-use algorithm
    if(alone_j){
      beta_aux[sample(N_aux,1),] <- beta[aux_s,] 
      if (!is.null(xi)) xi_aux[sample(seq_len(N_aux),1)] <- xi[aux_s] 
    }
    
    f_k <- double(N_aux + J)
    #The auxiliary are the easy ones
    fixed_aux <- tcrossprod(X[i, , ], beta_aux) ## matrix N x N_aux
    f_k[seq_len(N_aux)] <- colSums(dnorm(y[i, ], fixed_aux + w[i,], sqrt(sig2), 
                                            log = TRUE))
    

    fixed_expected <- tcrossprod(X[i, , ], beta) ## matrix N x J
    d_e <- dnorm(y[i, ], fixed_expected + w[i,], sqrt(sig2), log = TRUE)
    f_k[N_aux + 1:J] <-  if(J > 1) colSums(d_e) else sum(d_e)

    if (!is.null(xi)) {
      #The cluster labels require also the conditional densities of w
      Qw_1     <- drop(crossprod(w_demeaned[-i, ], Q[i, -i]))
      mu_c_aux <- tcrossprod(w_t_1[i, ], xi_aux) - s2_c[i] *  Qw_1
      mu_c     <- tcrossprod(w_t_1[i, ], xi)     - s2_c[i] *  Qw_1
      
      f_k[seq_len(N_aux)] <-  f_k[seq_len(N_aux)] + 
         colSums(dnorm(w[i,-1], mu_c_aux, sqrt(tau2*s2_c[i]), 
                                 log = TRUE))
      
      d_c <- dnorm(w[i,-1], mu_c, sqrt(tau2 * s2_c[i]),
            log = TRUE)

      f_k[N_aux + seq_len(J)]  <- f_k[N_aux + seq_len(J)] + 
        if (J > 1) colSums(d_c) else sum(d_c)
       
    }
    weight <- c(alpha * rep(1/N_aux, N_aux), nj)
    
    f_k <- exp(f_k-max(f_k)) * weight
    # print(f_k)
    f_k <- f_k/sum(f_k)
    
    hh <- sample.int(length(f_k), 1, prob = f_k)
    
    if (hh <= N_aux){
      #New dish at this table
      if (alone_j){
        #Same number of clusters
        nj[aux_s] <- 1
        beta[aux_s,] <- beta_aux[hh,]
        if (!is.null(xi))  xi[aux_s] <- xi_aux[hh]
      } else {
        J <- J + 1
        nj <- c(nj, 1)
        s[i] <- J
        beta <- rbind(beta, beta_aux[hh,])
        if (!is.null(xi)) xi <- c(xi, xi_aux[hh])
      }
      #Restore used auxiliary variable
      beta_aux[hh,] <-  mvnfast::rmvn(1, mu_beta, Sigma_beta)
      if (!is.null(xi)) xi_aux[hh] <- 2 * rbeta(1, a_xi, b_xi) - 1
    }else{
      #Old dish at this table
      h_k <- hh - N_aux
      nj[h_k] <- nj[h_k] + 1
      s[i] <- h_k
      if(alone_j){
        J <- J - 1
        s[s > aux_s] <- s[s > aux_s] - 1
        nj <- nj[-aux_s]
        beta <- beta[-aux_s,, drop=FALSE]
        if (!is.null(xi))  xi <- xi[-aux_s]
      }
    }
  }

  ###########################################################
  ## Update unique values for beta
  #if (is.null(dim(beta)[2])) beta <- matrix(beta, nrow = 1)
  beta.offset <- (y - offset - w - beta0)
  for(j in 1:J){
    index_j <- (1:I)[s == j]
    
    #Conjugate full-conditional for beta
    sumtXX <- matrix(0, p, p)
    aux_yX <- rep(0, p)
    for(i in index_j){
      sumtXX <- sumtXX + t(as.matrix(X[i, ,])) %*% as.matrix(X[i, , ])
      aux_yX <- aux_yX + beta.offset[i, ] %*% as.matrix(X[i, ,])
    }
    
    
    S_beta <- chol2inv(chol(Omega_beta + sumtXX/sig2))
    m_beta <- (mu_beta %*% Omega_beta + aux_yX/ sig2) %*% S_beta
    
    
    beta[j,] <- mvnfast::rmvn(1, mu = m_beta, sigma = S_beta)
  }
  # mu_all <- rep(mu_beta, J)  
  # Omega_all <- diag(rep(diag(Omega_beta), J))
  # s_alloc_vec <- rep(s, N)
  # Xs <- if (J > 1) model.matrix(~ - 1 + Xl:factor(s_alloc_vec)) else Xl
  # 
  # S_beta <- chol2inv(chol(Omega_all + crossprod(Xs)/sig2))
  # beta.offset2 <- crossprod(Xs, beta.offset)/sig2 + Omega_all %*% mu_all
  # m_beta      <- S_beta %*% beta.offset2
  # beta        <- drop(m_beta + crossprod(chol(S_beta), rnorm(length(m_beta))))
  
  # Sample unique values for xi
  if (!is.null(xi)) {
    for (j in seq_len(J)){
      index_j <- seq_len(I)[s == j]
      # MH step for xi
      #MH step for xi
      # xi_new <- log((1 + xi[j])/(1 - xi[j])) + sqrt(s_xi)*rnorm(1)
      # xi_new <- (exp(xi_new) - 1)/(1 + exp(xi_new))
      # 
      # #log-Acceptance rate
      # #Prior
      # log_accept <- a_xi * (log(1 + xi_new) - log(1 + xi[j])) + b_xi * (log(1-xi_new) - log(1-xi[j]))
      # #Likelihood
      # S2_c <- solve(Q[index_j,index_j])
      # # Mu_c <- -S2_c%*%Q[index_j,-index_j]%*%w[-index_j,1]
      # # log_accept <- log_accept + dmvn(w[index_j,1], Mu_c, tau2*S2_c, log = TRUE)
      # for(t in 2:N){
      #   Mu_c <- xi_new*w[index_j,t-1] - S2_c%*%Q[index_j,-index_j]%*%(w[-index_j,t] - xi[s[-index_j]]*w[-index_j,t-1])
      #   log_accept <- log_accept + dmvn(w[index_j,t], Mu_c, tau2*S2_c, log = TRUE)
      #   
      #   Mu_c <- xi[j]*w[index_j,t-1] - S2_c%*%Q[index_j,-index_j]%*%(w[-index_j,t] - xi[s[-index_j]]*w[-index_j,t-1])
      #   log_accept <- log_accept - dmvn(w[index_j,t], Mu_c, tau2*S2_c, log = TRUE)
      # }
      xi_new <- rtruncnorm(n = 1, a = -1, b = 1, mean = xi[j], sd = sqrt(s_xi))

      hastings <- 0 # = 0 because we truncate the normal symmetrically
                    # log(dtruncnorm(x = xi[j],  a = -1, b = 1, mean = xi_new,  sd = sqrt(s_xi))) -
                    # log(dtruncnorm(x = xi_new, a = -1, b = 1, mean = xi[j],   sd = sqrt(s_xi)))
      # log-Acceptance rate
      # Prior: scaled beta on [-1,1]: (x + 1)^a_xi * (1 - x)^b_xi
      log_accept <- hastings + a_xi * (log(1 + xi_new) - log(1 + xi[j])) +
        b_xi * (log(1 - xi_new) - log(1 - xi[j]))

      # Likelihood
      S2_c <- chol2inv(chol(Q[index_j, index_j]))

      SQW <- S2_c %*% Q[index_j,-index_j] %*% (w[-index_j, -1] - xi[s[-index_j]] * w_t_1[-index_j,])
      Mu_c1 <- xi_new * w_t_1[index_j, ] -  SQW
      Mu_c2 <- xi[j]  * w_t_1[index_j, ] -  SQW
      Mu_c1 <- as.matrix(Mu_c1)
      Mu_c2 <- as.matrix(Mu_c2)

      zeroIndex <- rep.int(0, length(index_j))
      cholS2_c <- chol(tau2 * S2_c)
      log_accept <- log_accept +
        sum(mvnfast::dmvn(t(w[index_j, -1] - Mu_c1), zeroIndex, cholS2_c, isChol = TRUE, log = TRUE)) -
        sum(mvnfast::dmvn(t(w[index_j, -1] - Mu_c2), zeroIndex, cholS2_c, isChol = TRUE, log = TRUE))

      accept <- 1
      if(is.nan(log_accept)){
        accept <- 0
      }else if( log_accept < 0 ){
        accept <- exp(log_accept)
      }

      if( runif(1) < accept ){
        xi[j] <- xi_new
      }
    }
  }

  return(list(s = s, beta = c(t(beta)), xi = xi))
}
