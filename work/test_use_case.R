library(Matrix)
library(mvtnorm)
library(slam)
library(spdep)
library(CARBayesSTBNP)
library(ggplot2)
##########################################
######## FIXING THE PARAMETERS ##########
##########################################
n1 <- 10L
n2 <- 10L


N <-  6 * 24 ## no of time points
rho <- 0.9
tau2 <- 1
P <- 2 # no of covariates
K <- 3
J <- 5L ## Number of clusters
sig2 <- 1 ## variance of the observations

######################
#### Simulation ######
######################
source("work/simulate_function.R")
sampleBeta0BNP <- FALSE
sim_obj <- simulate_spatio_temporal_model(n1 = n1, n2 = n2, N = N, J = J, K = K,
                                          rho = rho, tau2 = tau2, sig2 = sig2, 
                                          sampleBeta0BNP = sampleBeta0BNP,
                                          harmonic = FALSE)

df <- sim_obj$data
params <- sim_obj$params

formula <- sprintf("y ~ 1 + %s", paste0(colnames(df)[grepl("X", colnames(df))], collapse=" + "))

W <- params$W


niter  <- 200#00
burnin <- 50#000
thin   <- 1


####################
### CARar model
####################
library(mvnfast)
m_car_ar_bnp <- CARBayesSTBNP::ST.CARar(formula, data = df, W = W, family = "gaussian",
                                    n.sample = niter, burnin = burnin, thin = thin, 
                                    sampleBeta = "BNP", 
                                    sampleBNPgamma = TRUE,
                                    sampleBeta0BNP = FALSE
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

clobj <- mcclust::minbinder(mcclust::comp.psm(m_car_ar_bnp$samples$s))
table(clobj$cl, params$s)
