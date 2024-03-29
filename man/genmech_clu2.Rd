% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gendata.R
\name{genmech_clu2}
\alias{genmech_clu2}
\title{genmech clustering 2}
\usage{
genmech_clu2(npred = 10)
}
\arguments{
\item{npred}{number of predictive covariates used to generate the outcome}
}
\value{
a list of 8 elements.
\itemize{
\item \code{Y}: \eqn{n\times K \times}\code{nset} array. It contains \code{nset} matrices that store the outcome in the form of a multinomial experiment
\item \code{Yord}: \eqn{n\times}\code{nset} matrix. It contains the ordinal outcome for each replicated dataset
\item \code{treatment}: \eqn{n-}dimensional vector. It contains the treatment assigned
\item \code{pred}: \eqn{n\times}\eqn{(}\code{npred}\eqn{+}\code{nnoise}\eqn{)} matrix. It contains the predictive biomarkers (only the first \code{npred} are effectively used to generate the response)
\item \code{pred}: \eqn{n\times 2} matrix. It contains the prognostic biomarkers
\item \code{clu1}: contains the cluster label for patients assigned to Treatment 1
\item \code{clu2}: contains the cluster label for patients assigned to Treatment 2
\item \code{prob}: List of \eqn{T}. Each element is a \eqn{n\times K} matrix containing the response probabilities
}
}
\description{
Generates the data for the simulations scenarios S3 reported in the Supplementary Material.
It follows the strategy designed by Ma et al. (2019).
In particular, the outcome is a categorical variable representing \eqn{K=3} benefit-increasing levels.
Patients (\eqn{n}) are assigned to \eqn{T=2} competing treatments.
Moreover, the clustering is known and fixed. We generate predictive biomarkers as Scenarion 4 in Page and Quintana (2018)
}
\references{
Ma, J., Stingo, F. C., & Hobbs, B. P. (2019). Bayesian personalized
treatment selection strategies that integrate predictive with prognostic
determinants. \emph{Biometrical Journal}, \strong{61}(4), 902-917.
\url{https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700323}

Page, G. L. and Quintana, F. A. (2018). Calibrating covariate informed
product partition models. \emph{Statistics and Computing}, \strong{28}(5), 1009–1031.
}
