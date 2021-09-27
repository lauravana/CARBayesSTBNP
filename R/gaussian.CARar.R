gaussian.CARar <- function(formula, data=NULL, W, W.islands=NULL, burnin, n.sample, thin=1,  
                           prior.mean.beta=NULL, prior.var.beta=NULL, 
                           prior.nu2=NULL, prior.tau2=NULL, rho.S=NULL, rho.T=NULL, verbose=TRUE, 
                           sampleBeta=c("Normal", "BNP", "SFM"),
                           sampleBNPgamma = FALSE, 
                           sampleBeta0BNP = FALSE,
                           n.clusters = NULL, 
                           alpha_BNP = NULL,  N_aux = NULL) {
  ##############################################
  #### Format the arguments and check for errors
  ##############################################
  #### Verbose
  a <- common.verbose(verbose)  
  
  
  #### Frame object
  frame.results <- common.frame(formula, data, "gaussian")
  N.all <- frame.results$n
  X <- frame.results$X
  if (!sampleBeta0BNP) X <- X[,-1]
  p <- ncol(X)
  # X.sd <- frame.results$X.sd
  # X.mean <- frame.results$X.mean
  # X.indicator <- frame.results$X.indicator 
  offset <- frame.results$offset
  Y <- frame.results$Y
  which.miss <- frame.results$which.miss
  n.miss <- frame.results$n.miss  
  Y.DA <- Y      
  
  #### CAR quantities
  W.quants <- common.Wcheckformat.leroux(W)
  K <- I <-  W.quants$n
  N <- N.all / K
  W <- W.quants$W
  W.triplet <- W.quants$W.triplet
  W.n.triplet <- W.quants$n.triplet
  W.triplet.sum <- W.quants$W.triplet.sum
  n.neighbours <- W.quants$n.neighbours 
  W.begfin <- W.quants$W.begfin
  
  
  ### BNP

  if (is.null(n.clusters)) n.clusters <- min(25, round(K/2))
  J <- n.clusters
  if (is.null(alpha_BNP)) alpha_BNP <- 1
  if (is.null(N_aux)) N_aux <- 50
  
  #### Check on the rho arguments
  if(is.null(rho.S))
  {
    rho <- runif(1, 0, 0.9)
    fix.rho.S <- FALSE   
  }else
  {
    rho <- rho.S
    fix.rho.S <- TRUE
  }
  if (fix.rho.S) {
    if(!is.numeric(rho)) stop("rho.S is fixed but is not numeric.", call.=FALSE)  
    if(rho<0 ) stop("rho.S is outside the range [0, 1].", call.=FALSE)  
    if(rho>1 ) stop("rho.S is outside the range [0, 1].", call.=FALSE)    
  }
  if(is.null(rho.T))
  {
    if (!sampleBNPgamma) gamma <- runif(1) else gamma <- runif(n.clusters) 
    fix.rho.T <- FALSE   
  }else
  {
    gamma <- rho.T
    fix.rho.T <- TRUE
  }
  if (fix.rho.T) {
    if(!is.numeric(gamma)) stop("rho.T is fixed but is not numeric.", call.=FALSE)  
    if(gamma < 0) stop("rho.T is outside the range [0, 1].", call.=FALSE)  
    if(gamma > 1) stop("rho.T is outside the range [0, 1].", call.=FALSE)  
  }
  
  
  ## make an array
  X.array <- X
  dim(X.array) <- c(K, N, p)
  
  #### Priors
  if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
  if(is.null(prior.var.beta))  prior.var.beta  <- rep(1, p)
  if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
  if(is.null(prior.nu2))  prior.nu2 <- c(1, 0.01)
  prior.beta.check(prior.mean.beta, prior.var.beta, p)
  prior.var.check(prior.tau2)
  prior.var.check(prior.nu2)
  
  
  
  #### MCMC quantities - burnin, n.sample, thin
  common.burnin.nsample.thin.check(burnin, n.sample, thin)
  
  
  
  #############################
  #### Initial parameter values
  #############################
  if (!sampleBeta0BNP) {
    mod.glm <- glm(Y ~ X, offset = offset)
    beta0     <- mod.glm$coefficients[1]
    beta.mean <- mod.glm$coefficients[-1]   
    beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))[-1]
  } else {
    mod.glm <- glm(Y ~ 0 + X, offset = offset)
    beta0     <- 0
    beta.mean <- mod.glm$coefficients 
    beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
  }
  beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
  if (sampleBeta == "BNP" | sampleBNPgamma) {
    s_alloc <- sample.int(n.clusters, K, replace = TRUE)
    while (length(unique(s_alloc))< n.clusters) s_alloc <- sample.int(n.clusters, K, replace = TRUE)
    beta.mat <- matrix(rnorm(n.clusters * length(beta)), nrow = n.clusters, byrow = TRUE)
  }
  res.temp <- Y - beta0 - X %*% beta - offset
  res.sd <- sd(res.temp, na.rm=TRUE)
  phi <- rnorm(n=N.all, mean=0, sd = res.sd/2)
  tau2 <- var(phi)/2
  nu2 <- res.sd/2
  
  
  #### Matrix versions of quantities
  offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
  regression.mat <- matrix(X %*% beta, nrow=K, ncol=N, byrow=FALSE)
  phi.mat <- matrix(phi, nrow=K, ncol=N)   
  fitted <- as.numeric(beta0 + offset.mat + regression.mat + phi.mat)
  
  
  ###############################    
  #### Set up the MCMC quantities    
  ###############################
  #### Matrices to store samples
  n.keep <- floor((n.sample - burnin)/thin)
  samples.beta <-  vector("list", n.keep)
  samples.beta0 <-  array(NA, c(n.keep, 1)) 
  samples.s   <- array(NA, c(n.keep, I))
  samples.phi <- array(NA, c(n.keep, N.all))
  samples.J <- array(NA, c(n.keep, 1))
  samples.tau2 <- array(NA, c(n.keep, 1))
  samples.nu2 <- array(NA, c(n.keep, 1))
  if(!fix.rho.S) samples.rho <- array(NA, c(n.keep, 1))
  if(!fix.rho.T) samples.gamma <- array(NA, c(n.keep, ifelse(sampleBNPgamma, K, 1)))
  samples.fitted <- array(NA, c(n.keep, N.all))
  samples.loglike <- array(NA, c(n.keep, N.all))
  if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
  
  
  #### Specify the Metropolis quantities
  accept.all.rho <- accept.all.gamma <- rep(0,2)
  accept.rho <- accept.all.rho
  accept.gamma <- accept.all.gamma
  proposal.sd.rho <- 0.01
  proposal.sd.gamma <- 0.01
  tau2.shape <- prior.tau2[1] + N.all/2
  nu2.shape <- prior.nu2[1] + N.all/2        
  
  
  
  #############################
  #### Specify spatial elements
  #############################
  #### Spatial determinant
  if(!fix.rho.S) {
    Wstar <- diag(apply(W,1,sum)) - W
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.Q.W <-  0.5 * sum(log((rho * Wstar.val + (1-rho))))     
  }
  
  #### Beta update quantities
  data.precision.beta <- crossprod(X)
  if(length(prior.var.beta)==1)
  {
    prior.precision.beta <- 1 / prior.var.beta
  }else
  {
    prior.precision.beta <- diag(1/prior.var.beta)
  }
  
  
  #### Check for islands
  # if (is.null(W.islands)){
  #   W.list<- mat2listw(W)
  #   W.nb <- W.list$neighbours
  #   W.islands <- n.comp.nb(W.nb)
  # }
  # islands <- W.islands$comp.id
  # n.islands <- max(W.islands$nc)
  if(fix.rho.T) {
    if(rho==1 & gamma==1) {
      tau2.phi.shape <- prior.tau2[1] + prior.tau2[1] + ((N-1) * (K-1))/2
    }else if(rho==1)
    {
      tau2.phi.shape <- prior.tau2[1] + prior.tau2[1] + (N * (K-1))/2        
    }else if(gamma==1)
    {
      tau2.phi.shape <- prior.tau2[1] + prior.tau2[1] + ((N-1) * K)/2          
    }
  }
  
  ## needed for BNP part
  hyperparameters_P0 <- list(mu_beta = prior.mean.beta, 
                             Omega_beta = prior.precision.beta,
                             Sigma_beta = diag(prior.var.beta),
                             a_xi = 1, 
                             b_xi = 1,
                             s_xi =proposal.sd.gamma^2)
  
  ###########################
  #### Run the Bayesian model
  ###########################
  #### Start timer
  if(verbose) {
    cat("Generating", n.keep, "post burnin and thinned (if requested) samples.\n", sep = " ")
    progressBar <- txtProgressBar(style = 3)
    percentage.points<-round((1:100/100)*n.sample)
  } else {
    percentage.points<-round((1:100/100)*n.sample)     
  }
  
  #j<-1
  #### Create the MCMC samples     
  for(j in 1:n.sample) {
    
    ####################################
    ## Sample from Y - data augmentation
    ####################################
    if (n.miss > 0) {
      Y.DA[which.miss==0] <- rnorm(n = n.miss, mean = fitted[which.miss==0], sd = sqrt(nu2))    
    }
    Y.DA.mat <- matrix(Y.DA, nrow=I, ncol = N, byrow=FALSE)
    
    
    
    ##################
    ## Sample from nu2
    ##################
    nu2.offset <- as.numeric(Y.DA.mat - beta0 - regression.mat - offset.mat -  phi.mat)
    nu2.scale  <- prior.nu2[2]  + sum(nu2.offset ^ 2)/2
    nu2        <- 1 / rgamma(1, nu2.shape, nu2.scale)

    ############################
    ## Sample from beta & xi ##
    ###########################
    
    if (sampleBeta == "Normal") {
      fc.precision <- prior.precision.beta + data.precision.beta / nu2
      fc.var <- solve(fc.precision)
      beta.offset <- as.numeric(Y.DA.mat - beta0 - offset.mat - phi.mat)
      beta.offset2 <- crossprod(X, beta.offset) / nu2 + prior.precision.beta %*% prior.mean.beta
      fc.mean <- fc.var %*% beta.offset2
      chol.var <- chol(fc.var)
      beta <- drop(fc.mean + crossprod(chol.var, rnorm(p)))
      regression.mat <- matrix(X %*% beta, nrow=K, ncol=N, byrow=FALSE)
    }
    if (sampleBeta == "BNP") {
      Sys.time()
      out_BNP <- sample_BNP(y = Y.DA.mat, Xl = X, X = X.array, 
                            w = phi.mat, offset = offset.mat, beta0=beta0,  s = s_alloc, 
                            beta = beta.mat, sig2 = nu2, tau2=tau2, 
                            xi = if (sampleBNPgamma) gamma else NULL,
                            rho = rho,
                            alpha = alpha_BNP, N_aux=N_aux,
                            hyperparameters_P0 = hyperparameters_P0)
      Sys.time()
      betavec <- out_BNP$beta
      s_alloc <-  out_BNP$s
      J <- length(unique(s_alloc))
      beta.mat <- matrix(betavec, nrow = max(s_alloc), byrow = TRUE)
      s_alloc_vec <- rep(s_alloc, N)
      Xs <- if (J > 1)  model.matrix(~ - 1 + X:factor(s_alloc_vec)) else X
      regression.mat <- matrix(Xs %*% betavec, nrow=K, ncol=N, byrow=FALSE)  
    }
    
    if (sampleBeta == "BNP" & sampleBNPgamma) {
      gamma <- out_BNP$xi
      gamma.vec <- gamma[s_alloc]
    } else {
      gamma.vec <- rep(gamma, K)
    }
    
    ##########################
    ## Sample from beta0 ##
    ###########################
    if (sampleBeta0BNP) {
      beta0 <- 0
    } else {
      fc0.precision <- 0.01 + N * K/ nu2
      fc0.var <- 1/fc0.precision
      beta0.offset <- as.numeric(Y.DA.mat - regression.mat - offset.mat - phi.mat)
      
      beta0.offset2 <-  sum(beta0.offset)/nu2
      fc0.mean <- fc0.var * beta0.offset2
      beta0 <- fc0.mean + sqrt(fc0.var) * rnorm(1)
    }

  
    # 
    ####################
    ## Sample from phi
    ####################
    phi.offset <- Y.DA.mat - offset.mat - regression.mat - beta0
    den.offset <- rho * W.triplet.sum + 1 - rho
    phi.temp <- gaussianarcarupdate(W.triplet, W.begfin, W.triplet.sum,  K, N, phi.mat, tau2, nu2, gamma.vec, rho, phi.offset, den.offset)
    phi <- as.numeric(phi.temp) - mean(as.numeric(phi.temp))
    phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)



    ####################
    ## Sample from gamma
    ####################
    if(!fix.rho.T & !sampleBNPgamma) {
      temp2 <- gammaquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, rho)
      mean.gamma <- temp2[[1]] / temp2[[2]]
      sd.gamma <- sqrt(tau2 / temp2[[2]])
      gamma <- rtruncnorm(n=1, a=-1, b=1, mean=mean.gamma, sd=sd.gamma)
      gamma.vec <- rep(gamma, K)
    }

    
    ####################
    ## Samples from tau2
    ####################
    
    temp3 <- tauquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, rho, gamma.vec)
    tau2.scale <- temp3 + prior.tau2[2]
    tau2 <- 1 / rgamma(1, tau2.shape, scale=(1/tau2.scale)) 

    ##################
    ## Sample from rho
    ##################
    if(!fix.rho.S)
    {
      proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)
      temp4 <- tauquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, proposal.rho, gamma.vec)
      det.Q.W.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))
      logprob.current <- N * det.Q.W - temp3 / tau2
      logprob.proposal <- N * det.Q.W.proposal - temp4 / tau2
      hastings <- log(dtruncnorm(x=rho, a=0, b=1, mean=proposal.rho, sd=proposal.sd.rho)) -
        log(dtruncnorm(x=proposal.rho, a=0, b=1, mean=rho, sd=proposal.sd.rho)) 
      prob <- exp(logprob.proposal - logprob.current + hastings)
      if(prob > runif(1))
      {
        rho <- proposal.rho
        det.Q.W <- det.Q.W.proposal
        accept.rho[1] <- accept.rho[1] + 1           
      }else
      {}              
      accept.rho[2] <- accept.rho[2] + 1       
    }else
    {}
    


    #########################
    ## Calculate the deviance
    #########################
    fitted <- as.numeric(offset.mat + regression.mat + beta0 + phi.mat)
    loglike <- dnorm(Y, mean = fitted, sd = rep(sqrt(nu2), N.all), 
                     log=TRUE)
    
    
    
    ###################
    ## Save the results
    ###################
    if(j > burnin & (j-burnin)%%thin==0) {
      ele <- (j - burnin) / thin
      samples.beta[[ele]] <- switch(sampleBeta,
                                     "Normal" = beta, 
                                     "BNP"    = beta.mat) 
      if (sampleBeta != "Normal"|sampleBNPgamma) samples.s[ele, ] <- s_alloc
      if (sampleBeta != "Normal"|sampleBNPgamma) {
        samples.J[ele, ] <- length(unique(s_alloc))
      }
      samples.phi[ele, ] <- as.numeric(c(phi.mat))

      if(!fix.rho.S) samples.rho[ele, ] <- rho
      if(!fix.rho.T) samples.gamma[ele, ] <- if (sampleBNPgamma) gamma.vec else gamma
      samples.tau2[ele, ] <- tau2
      samples.nu2[ele, ] <- nu2
      samples.beta0[ele, ] <- beta0
      samples.fitted[ele, ] <- fitted
      samples.loglike[ele, ] <- loglike
      if(n.miss>0) samples.Y[ele, ] <- Y.DA[which.miss==0]
    }else
    {}
    
    
    
    #########################################
    ## Self tune the acceptance probabilities
    #########################################
    k <- j/100
    if(ceiling(k)==floor(k))
    {
      if(!fix.rho.S) proposal.sd.rho <- common.accceptrates2(accept.rho[1:2], proposal.sd.rho, 40, 50, 0.5)
     # if(!fix.rho.T & !sampleBNPgamma) proposal.sd.gamma <- common.accceptrates2(accept.gamma[1:2], proposal.sd.gamma, 40, 50, 0.5)
      accept.all.rho <- accept.all.rho + accept.rho
      accept.rho <- rep(0,2)
      #accept.all.gamma <- accept.all.gamma + accept.gamma
      #accept.gamma <- rep(0,2)
    }else
    {}
    
    
    
    ################################       
    ## print progress to the console
    ################################
    if(j %in% percentage.points & verbose)
    {
      setTxtProgressBar(progressBar, j/n.sample)
    }
  }
  
  
  #### end timer
  if(verbose)
  {
    cat("\nSummarising results.")
    close(progressBar)
  }else
  {}
  
  
  
  ###################################
  #### Summarise and save the results 
  ###################################
  #### Compute the acceptance rates
  if(!fix.rho.S)
  {
    accept.rho <- 100 * accept.all.rho[1] / accept.all.rho[2]
  }else
  {
    accept.rho <- NA    
  }
  if(!fix.rho.T & !sampleBNPgamma)
  {
    accept.gamma <- 100 #* accept.all.gamma[1] / accept.all.gamma[2]
  }else
  {
    accept.gamma <- NA    
  }
  accept.phi  <- 100
  accept.final <- c(accept.phi, accept.rho, accept.gamma)
  names(accept.final) <- c("phi", "rho.S", "rho.T")
  
  
  #### Compute the fitted deviance
  ## 1) Deal with betas
  if (sampleBeta == "Normal") {
    mean.beta <- colMeans(do.call("rbind", samples.beta))
    regression.mat <- matrix(X %*% mean.beta, nrow=K, ncol=N, byrow=FALSE)   
  }
  if (sampleBeta == "BNP") {

    uniq <- unique(samples.J)
    K0 <- uniq[which.max(tabulate(match(samples.J, uniq)))[1]]
    print(K0)
    # clobj <- mcclust::minbinder(mcclust::comp.psm(samples.s[samples.J == K0,]))
    # s_post <- clobj$cl # posterior allocation
    s_post <- apply(samples.s[samples.J == K0,],2,function(x) {
      uniq <- unique(x)
      uniq[which.max(tabulate(match(x, uniq)))[1]]
    })
    
    samples.beta.tmp <- samples.beta[sapply(samples.beta, function(x) nrow(x) == K0)]
    mean.beta <- Reduce("+", samples.beta.tmp)/length(samples.beta.tmp)
    beta.tmp <- if (K0 == 1) drop(mean.beta) else c(t(mean.beta))
    Xs <- if (K0 > 1) model.matrix(~ - 1 + X : factor(rep(s_post, N))) else X
    regression.mat <- matrix(Xs %*% beta.tmp, nrow = I, ncol = N, byrow=FALSE)   
  }
  mean.phi   <- matrix(apply(samples.phi, 2, mean), nrow=K, ncol=N)
  mean.beta0 <- mean(samples.beta0)
  fitted.mean <- as.numeric(offset.mat + mean.phi + regression.mat + mean.beta0)
  nu2.mean <- mean(samples.nu2)
  deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.mean, sd = rep(sqrt(nu2.mean),N.all), log = TRUE), na.rm=TRUE)

  #### Model fit criteria
  modelfit <- common.modelfit(samples.loglike, deviance.fitted)
   
  
  #### Create the fitted values and residuals
  # fitted.values <- apply(samples.fitted, 2, mean)
  # response.residuals <- as.numeric(Y) - fitted.values
  # pearson.residuals <- response.residuals /sqrt(nu2.mean)
  # residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)
  # 
  # 
  #### Transform the parameters back to the origianl covariate scale.
  # samples.beta.orig <- common.betatransform(samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
  # 
  # if (sampleBNP) {
  #   #samples.beta.all <- do.call("cbind", lapply(1:length(samples.beta), function(i) 
  #   #  samples.beta[[i]][samples.s[[i]]]))
  #   ## todo transformation
  #   #samples.beta.orig <- t(samples.beta.all)
  # } else {
  #   samples.beta.all <- t(do.call("cbind", samples.beta))
  #   samples.beta.orig <- common.betatransform(samples.beta.all, X.indicator, X.mean, X.sd, p, FALSE)
  #   samples.beta.orig <- mcmc(samples.beta.orig)
  # }
  #### Create a summary object
  # 
  # summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
  # summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(100,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
  # rownames(summary.beta) <- colnames(X)
  # colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
  # 
  # summary.hyper <- array(NA, c(4, 7))     
  # summary.hyper[1,1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
  # summary.hyper[2,1:3] <- quantile(samples.nu2, c(0.5, 0.025, 0.975))
  # rownames(summary.hyper) <- c("tau2", "nu2", "rho.S", "rho.T")     
  # summary.hyper[1, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.tau2)), geweke.diag(mcmc(samples.tau2))$z)     
  # summary.hyper[2, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.nu2)), geweke.diag(mcmc(samples.nu2))$z)     
  # 
  # if(!fix.rho.S)
  # {
  #   summary.hyper[3,1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
  #   summary.hyper[3, 4:7] <- c(n.keep, accept.rho, effectiveSize(mcmc(samples.rho)), geweke.diag(mcmc(samples.rho))$z)  
  # }else
  # {
  #   summary.hyper[3, 1:3] <- c(rho, rho, rho)
  #   summary.hyper[3, 4:7] <- rep(NA, 4)
  # }
  # if(!fix.rho.T)
  # {
  #   summary.hyper[4,1:3] <- quantile(samples.gamma, c(0.5, 0.025, 0.975))
  #   summary.hyper[4, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.gamma)), geweke.diag(mcmc(samples.gamma))$z)  
  # }else
  # {
  #   summary.hyper[4, 1:3] <- c(gamma, gamma, gamma)
  #   summary.hyper[4, 4:7] <- rep(NA, 4)
  # }   
  # 
  # summary.results <- rbind(summary.beta, summary.hyper)
  # summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
  # summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
  # 
  # 
  #### Compile and return the results
  #### Harmonise samples in case of them not being generated
  if (sampleBNPgamma) {
    if(!fix.rho.S) {
      samples.rhoext <- samples.rho
      names(samples.rhoext) <- "rho.S"
    } else {
      samples.rhoext <- NA
    }
  } else {
    if(!fix.rho.S) {
      if (!fix.rho.T) {
        samples.rhoext <- cbind(samples.rho, samples.gamma)
        colnames(samples.rhoext) <- c("rho.S", "rho.T")
      } else {
        samples.rhoext <- samples.rho
        names(samples.rhoext) <- "rho.S"
      }
    } else {
      if (!fix.rho.T) {
        samples.rhoext <- samples.gamma
        names(samples.rhoext) <- "rho.T"
      } else {
        samples.rhoext <- NA
      }
    }
  }
  
  if(n.miss==0) samples.Y = NA
  
  samples <- list(beta=(samples.beta),
                  s = mcmc(samples.s), 
                  J =if (sampleBeta != "Normal") mcmc(samples.J) else NULL, 
                  phi=mcmc(samples.phi), 
                  rho=mcmc(samples.rhoext), 
                  gamma=if (sampleBNPgamma) mcmc(samples.gamma) else NULL, 
                  tau2=mcmc(samples.tau2), 
                  nu2=mcmc(samples.nu2), 
                  beta0=mcmc(samples.beta0),
                  fitted=mcmc(samples.fitted), 
                  Y=mcmc(samples.Y))
  
  model.string <- c("Likelihood model - Gaussian (identity link function)", "\nLatent structure model - Autoregressive CAR model\n")
  results <- list(#summary.results=summary.results, 
    samples=samples, 
    #fitted.values=fitted.values, 
    #residuals=residuals,
    modelfit=modelfit, 
    accept=accept.final, #localised.structure=NULL, 
    formula=formula, model=model.string)
  class(results) <- "CARBayesST"
  
  #### Finish by stating the time taken 
  if(verbose)
  {
    b<-proc.time()
    cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
  }else
  {}
  return(results)
}
