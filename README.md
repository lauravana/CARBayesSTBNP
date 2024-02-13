# CARBayesSTBNP R package

This R package provides functionality for estimating a CAR-AR linear model with a 
Bayesian non-parametric (BNP) prior on the area-specific regression coefficients. 
This results in a clustering of the areal units based on their regression 
coefficients.

You can install the package using `remotes::install_github()`:

```
if (!require("remotes")) install.packages("remotes")
library("remotes")
install_github("lauravana/CARBayesSTBNP")
```

