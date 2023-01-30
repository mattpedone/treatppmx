
<!-- README.md is generated from README.Rmd. Please edit that file -->

# treatppmx

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/treatppmx)](https://CRAN.R-project.org/package=treatppmx)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

An R package for personalized treatment selection at a single decision
point.

Bayesian predictive model for personalized treatment selection for new
untreated patients, which leverages known predictive and prognostic
biomarkers. In particular, predictive biomarkers are exploited to inform
a product partition model with covariates (PPMx) to obtain homogeneous
clusters.

The implementation has been done in C++ through the use of Rcpp and
RcppArmadillo.

Authors: Matteo Pedone, Raffaele Argiento, Francesco Stingo

Maintainer: Matteo Pedone.

## Installation

You can install the development version of treatppmx from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mattpedone/treatppmx")
```

**NOTE** that this package depends on
[vweights](https://github.com/mattpedone/vweights) and
[mcclust.ext](https://github.com/sarawade/mcclust.ext).

It has only been tested on a PC running Ubuntu 20.04.2 LTS.

## References 

We describe the method implemented in treatppmx in a preprint available at [arxive](https://arxiv.org/pdf/2210.06030.pdf).
