# PMLE4SCR

This R package implements the two-stage pseudo maximum likelihood estimation (PMLE) for the copula-based analysis of the semi-competing risks data. The marginal distributions are modeled by semiparametric transformation regression models, and the dependence between the bivariate event times is specified by a parametric copula function with the copula parameter that may depend on some covariates.

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
## Example

```{r}
data(BMT, package = "SemiCompRisks")
data = BMT %>%
        mutate(g = factor(g, levels = c(2, 3, 1),
                    labels = c("AML-low", "AML-high", "ALL")))
myfit = PMLE4SCR(data, time = "T2", death = "T1",
                        status_time = "delta2", status_death = "delta1",
                        T.fmla = ~ g, D.fmla = ~ g,
                        copula.family = "Clayton",
                        copula.control = list(link = "identity", formula = ~ g),
                        initial = c(2, 0, 0))
myfit$PMLE$gamma
myfit$PMLE$betaT
```
