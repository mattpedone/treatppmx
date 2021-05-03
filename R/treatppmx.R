#' Treatment PPMx
#'
#' The implementation has been done in \code{C++} through the use of \code{Rcpp} and \code{RcppArmadillo}.
#' @author
#' matt
#'
#' @docType package
#' @name treatppmx
#'
#' @useDynLib treatppmx, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom mcclust comp.psm
#' @importFrom mcclust.ext minbinder.ext
#' @importFrom mcclust.ext minVI
#' @importFrom mclust adjustedRandIndex
#' @importFrom coda effectiveSize
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom stats acf
#' @importFrom graphics par
#' @importFrom mvtnorm rmvnorm
#' @importFrom multiROC multi_roc
#' @importFrom multiROC plot_roc_data
#' @importFrom pROC multiclass.roc
NULL
