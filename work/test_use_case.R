library(Matrix)
library(mvtnorm)
library(slam)
library(spdep)
library(CARBayesSTBNP)
library(ggplot2)
##########################################
######## FIXING THE PARAMETERS ##########
##########################################
n1 <- 3L
n2 <- 3L


N <-  6 * 24 ## no of time points
rho <- 0.9
tau2 <- 1
P <- 2 # no of covariates
K <- 3
J <- 2L ## Number of clusters
sig2 <- 1 ## variance of the observations

######################
#### Simulation ######
######################
source("work/simulate_function.R")
sim_obj <- simulate_spatio_temporal_model(n1 = n1, n2 = n2, N = N, J = J, K = K,
                                          rho = rho, tau2 = tau2, sig2 = sig2, 
                                          xi = c(0.5,0.9),
                                          harmonic = FALSE, 
                                          #freq = freq, dummy_phit_lev = 1, 
                                          hyperparameters_P0 = list(mu_beta = rep(0, K),
                                                                    Sigma_beta = diag(K),
                                                                    alpha_xi = 0.5,
                                                                    beta_xi  = 0.5))

df <- sim_obj$data
params <- sim_obj$params


formula <- sprintf("y ~ 1+ %s", paste0(colnames(df)[grepl("X", colnames(df))], collapse=" + "))
W <- as.matrix(read.table(sprintf("adj_matrix_n1_%i_n2_%i_T_%i.csv", 
                                  n1, n2, N)))

#data = df; W=W; family = "gaussian";n.sample = 10; burnin = 1; thin= 1;

niter  <- 100#00
burnin <- 50#000
thin   <- 1


####################
### CARar model
####################
library(mvnfast)
m_car_ar_bnp <- CARBayesSTBNP::ST.CARar(formula, data = df, W=W, family = "gaussian",
                                    n.sample = niter, burnin = burnin, thin = thin, 
                                    sampleBeta = "BNP", sampleBNPgamma = TRUE,
                                    sampleBeta0BNP = TRUE
)

head(m_car_ar_bnp$samples$beta0)

head(m_car_ar_bnp$samples$beta)

head(m_car_ar_bnp$samples$rho)

head(m_car_ar_bnp$samples$rho)
mean(m_car_ar_bnp$samples$rho[,1])
mean(m_car_ar_bnp$samples$rho[,2])
plot(m_car_ar_bnp$samples$rho[-c(1:1000),1])#,ylim=c(0.7,1))
plot(m_car_ar_bnp$samples$rho[-c(1:1000),2])
plot(m_car_ar_bnp$samples$nu2)
plot(m_car_ar_bnp$samples$tau2)
mean(m_car_ar_bnp$samples$tau2)
mean(m_car_ar_bnp$samples$nu2[-c(1:100)])