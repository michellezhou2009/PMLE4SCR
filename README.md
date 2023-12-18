# PMLE4SCR

This R package implements the two-stage pseudo maximum likelihood estimation (PMLE) for the copula-based analysis of the semi-competing risks data. The marginal distributions are modeled by semiparametric transformation regression models, and the dependence between the bivariate event times is specified by a parametric copula function with the copula parameter that may depend on some covariates

## Instrall `PMLE4SCR`

To install this R package, you need to first install the `devtool` package via
```{r}
install.packages("devtools")
```
To install the `PMLE4SCR` package,
```{r}
library(devtools)
install_github("michellezhou2009/PMLE4SCR")
library(PMLE4SCR)
```
