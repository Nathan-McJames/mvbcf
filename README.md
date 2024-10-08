
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mvbcf <img src="man/figures/logo.png" align="right" height="120" alt="" />

<!-- badges: start -->
<!-- badges: end -->

Welcome to the mvbcf R package! This R package implements the
Multivariate Bayesian Causal Forest model introduced in McJames et al.,
2024:

McJames, N., O’Shea, A., Goh, Y. C. & Parnell, A. (2024). Bayesian
causal forests for multivariate outcomes: Application to Irish data from
an international large scale education assessment. *Journal of the Royal
Statistical Society Series A: Statistics in Society*, 1-23.
<https://doi.org/10.1093/jrsssa/qnae049>.

See below for details on installation. For worked examples and
quick-start guides showing how to use the package, check out the
vignettes and other materials. For questions, suggestions, or other
queries, please reach out to me at mcjamesnathan@yahoo.ie.

## Installation

You can install the development version of mvbcf with the following:

``` r

if (!require("devtools")) {
  install.packages("devtools")
}
install_github("Nathan-McJames/mvbcf")
#Or to install vignettes with the package as well:
#install_github("Nathan-McJames/mvbcf", build_vignettes = TRUE)
```
