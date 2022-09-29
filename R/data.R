#' 152 simulated patients data.
#'
#' A data frame 92 genes coming from 152 patients. This data were obtained
#' to emulate the dependence structure arising in sequencing data.
#' The simulated scenario starts from a well-known dataset of leukemia. The
#' process to obtained these data from the original Leukemia dataset is extensively
#' detailed in Ma et al. (2019).
#'
#' @references
#' Ma, J., Stingo, F. C., & Hobbs, B. P. (2019). Bayesian personalized
#' treatment selection strategies that integrate predictive with prognostic
#' determinants. \emph{Biometrical Journal}, \strong{61}(4), 902-917.
#' \url{https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700323}
#'
#' The Leukemia dataset containing gene expression levels for a total of
#' 5,000 genes across 38 patients.
#'
#' @format A data frame containing 152 observation of 92 variables:
#' \describe{
#'   \item{1-90}{predictive covariates}
#'   \item{91-92}{prognostic covariates}
#'          }
#'
#' @source \url{https://www.pnas.org/content/101/12/4164.full?tab=ds}
#'
"simupats"
