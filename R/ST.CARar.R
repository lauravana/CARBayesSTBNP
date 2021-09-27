ST.CARar <- function(formula, family, data=NULL, 
                     W, W.islands=NULL, burnin, n.sample, thin=1, 
                     prior.mean.beta=NULL, prior.var.beta=NULL, 
                     prior.nu2=NULL, prior.tau2=NULL,
                     rho.S=NULL, rho.T=NULL, verbose=TRUE, MALA=FALSE,
                     sampleBeta=c("Normal", "BNP", "SFM"), sampleBNPgamma=FALSE,
                     sampleBeta0BNP = FALSE,
                     n.clusters = NULL, alpha_BNP = NULL,  N_aux = NULL
                    )
{
    ## This is a wrapper function
    if(is.null(family)) stop("the family argument is missing", call.=FALSE)
    #### Run the appropriate model according to the family argument
    if(family=="gaussian") {
        model <- gaussian.CARar(formula = formula, data = data,  
                                W = W, W.islands = W.islands, 
                                burnin = burnin, n.sample = n.sample, thin = thin,
                                prior.mean.beta=prior.mean.beta, 
                                prior.var.beta=prior.var.beta, 
                                prior.nu2 = prior.nu2, prior.tau2 = prior.tau2,
                                rho.S = rho.S, rho.T = rho.T, 
                                sampleBeta = sampleBeta, 
                                sampleBeta0BNP = sampleBeta0BNP,
                                sampleBNPgamma = sampleBNPgamma,
                                verbose = verbose,
                                n.clusters = n.clusters, alpha_BNP = alpha_BNP,  N_aux = N_aux)          
    } else {
        stop("The family argument is not `gaussian'.", call.=FALSE)     
    }
    return(model)     
}
