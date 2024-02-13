\name{CARBayesSTBNP-package}
\alias{CARBayesSTBNP-package}
\alias{CARBayesSTBNP}
\docType{package}
\title{
  Bayesian spatio-temporal model for areal data with normal and BNP priors on regression coefficients
}
\description{
  Bayesian spatio-temporal for areal data with normal and BNP priors on regression coefficients. Adapted from the CARBayesST package and from Cadonna et al (2019).
}
\details{
}
\author{
Maintainer: Laura Vana <laura.vana@tuwien.ac.at>.
}
\references{
  
}
\keyword{ package }
\seealso{

}
\examples{
  \dontrun{
     data("data_car_ar_bnp", package = "CARBayesSTBNP")
     W <- as.matrix(build_Q_sparse(10, 10, rho = 0.9)$W)
     m_car_ar_bnp_2 <- CARBayesSTBNP::ST.CARar(
       y ~ 1 + X1 + X2 + X3, 
       data = data_car_ar_bnp, W = W, 
       family = "gaussian",
       n.sample = 1000, 
       burnin = 500, 
       thin = 1, 
       sampleBeta = "BNP", sampleBNPgamma = FALSE, sampleBeta0BNP = FALSE)
  }
}
