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
#' @importFrom mclust adjustedRandIndex
#' @importFrom coda effectiveSize
# @importFrom mvtnorm rmvnorm
NULL