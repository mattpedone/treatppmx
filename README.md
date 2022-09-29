# treatppmx

An R package for personalized treatment selection at a single decision point.

Bayesian predictive model for personalized treatment selection for new
untreated patients, which leverages known predictive and prognostic biomarkers.
In particular, predictive biomarkers are exploited to inform a product partition model with covariates (PPMx) to obtain homogeneous clusters.

The implementation has been done in C++ through the use of Rcpp and RcppArmadillo.

Authors: Matteo Pedone, Raffaele Argiento, Francesco Stingo

Maintainer: Matteo Pedone.

## Installation

You can install the development version from [GitHub](https://CRAN.R-project.org), using the **devtools** package with:

``` r
devtools::install_github("mattpedone/treatppmx")
```

**NOTE** that this package depends on [vweights](https://github.com/mattpedone/vweights). 

It has only been tested on a PC running Ubuntu 20.04.2 LTS.

## References

TBA
