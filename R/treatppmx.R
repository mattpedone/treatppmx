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
# @importFrom multiROC multi_roc
# @importFrom multiROC plot_roc_data
#' @importFrom vweights computev
#' @importFrom stats cor.test
#' @importFrom stats cor
#' @importFrom stats kmeans
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 geom_raster
#' @importFrom ggplot2 scale_fill_continuous
#' @importFrom ggplot2 geom_tile
#' @importFrom reshape2 melt
NULL
