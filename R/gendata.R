##################################################################################
## Data generating mechanism for simulation scenarios
## It follows the code by Junsheng Ma. I just re-arranged it
## A total of 92 genes are selected for 152 simuluated patients               ####
## Also simulate 100 realizations data sets                                   ####
##################################################################################

#' genoutcome
#'
#' Generates outcome variables. The response is an ordinal outcome variable.
#'
#' See Ma et al. (2019) Supplementary Material for all the details.
#'
#'
#' @param nobs number of observations
#' @param alpha intercept
#' @param beta1 coefficient for the predictive covariates. More precisely it is
#' a coefficient for the combination of the pca to generate the exponetial.
#' @param beta2 coefficient for the first prognostic covariate (x2)
#' @param beta3 coefficient for the second prognostic covariate (x3)
#' @param metx combination of the pca to generate the exponetial
#' @param x2 first prognostic covariate (x2)
#' @param x3 second prognostic covariate (x3)
#' @return a \eqn{\code{nobs}\times 4} matrix. The first column contains the generated
#' outcome variable, the other three the probabilities for the categories.

genoutcome <- function(nobs, alpha, beta1, beta2, beta3, metx, x2, x3){
  prob <- matrix(0, nrow = nobs, ncol = 3)
  myy <- rep(0, nobs)

  for(i in 1:nobs){
    temp <- NULL
    temp <- beta1*metx[i] + beta2*x2[i] + beta3*x3[i]
    eta1 <- alpha[1] + temp[1]
    eta2 <- alpha[2] + temp[2]

    prob[i,1] <- 1/((exp(eta1)+1)*(exp(eta2)+1))
    prob[i,2] <- exp(eta1)/((exp(eta1)+1)*(exp(eta2)+1))
    prob[i,3] <- exp(eta2)/(exp(eta2)+1)
    myy[i] <- sample(c(0,1,2), size = 1, prob = prob[i,])
  }
  outcome <- round(cbind(myy,prob), 4)
  return(outcome)
}

#' genmech
#'
#' Generates the data for the different simulations scenarios. It follows the
#' strategy designed by Ma et al. (2019). It is used to generate the data for all
#' the simulation study scenarios we used. In particular, it has three benefit-incresing levels. Patients are assigned
#' to two competing treatments.
#' The ordinal outcome variable are generated using two separate
#' continuation-ratio logistic functions. The information carried by many
#' predictive features can be summarized via dimension reduction techniques.
#' This implementation, in particular, considers thr first two principal
#' components of a PCA.
#'
#' See Ma et al. (2019) Supplementary Material for all the details.
#'
#' @param npred number of predictive covariates used to generate the outcome
#' @param mydata input dataset file
#' @param progscen Prognostic covariates option:
#'   1 - Prognostic biomarkers are considered in the original scale (default)
#'   2 - Prognostic biomarkers are transformed
#' @param predscen Predictive covariates option:
#'   1 - 10 Predictive biomarkers are considered to generate outcomes (default)
#'   2 - \code{nnoise} noisy std normals are added to the design matrix of predictive
#'   covariates in addition to 10 the predictive biomarkers considered to
#'   generate outcomes
#' @param nnoise number of noisy covariates added to predictive biomarkers
#' @param nset number of replicated scenarios generated
#' @return a list of 5 elements.
#'   \itemize{
#'   \item 1 - List of \code{nset} outcome variables.
#'   \item 2 - List of \code{nset} outcome variables in ordinal notation.
#'   \item 3 - A vector of assigned treatments.
#'   \item 4 - A matrix of all biomarkers.
#'   \item 5 - A list of outcome probabilities (one element for each treatment).
#' }
#' @export

genmech <- function(npred = 10, mydata = "data/Simu152pats.txt", progscen = 1,
                    predscen = 1, nnoise = 15, nset = 30){

  genenorm <- t(scale(t(read.table(file = mydata))))
  mypca <- prcomp(t(genenorm[c(1:npred),]))

  nobs <- length(mypca$x[,1])

  # Predictive Markers
  ## use the combination of the pca to generate the exponential
  metx1 <- scale(mypca$x[,1]) + scale(mypca$x[,2]) + 0.50
  metx <- sign(metx1)*(sign(metx1)*metx1)^(0.2) * 0.45
  ## for the noisy scenario

  # Prognostic Markers
  ## original scale
  x2 <- genenorm[91,]
  x3 <- genenorm[92,]
  ## transformation
  myx2 <- sign(x2)*(sign(x2)*x2)^(0.5)
  myx3 <- sign(x3)*(sign(x3)*x3)^(0.2)

  # pmts probabilities for treatment 1
  alpha1 <- c(-0.5, -1)
  beta11 <- c(2, 2.6)
  # pmts probabilities for treatment 2
  alpha2 <- c(0.7, -1)
  beta21 <- c(-1, -3)
  # pmts probabilities with prognostic only
  alpha3 <- c(1, -0.5)
  beta2 <- c(1, 0.5)
  beta3 <- c(0.7,1)

  if(progscen == 1){
    #probabilities for treatment 1
    prob1 <- genoutcome(nobs, alpha1, beta11 ,c(0,0), c(0,0), metx, x2, x3)
    #probabilities for treatment 2
    prob2 <- genoutcome(nobs, alpha2, beta21, c(0,0), c(0,0), metx, x2, x3)
    probprog <- genoutcome(nobs, alpha3, c(0,0), beta2, beta3, metx, x2, x3)
    Zprogcov <- cbind(x2, x3)
  }
  if(progscen == 2){
    #probabilities for treatment 1
    prob1 <- genoutcome(nobs, alpha1, beta11 ,c(0,0), c(0,0), metx, myx2, myx3)
    #probabilities for treatment 2
    prob2 <- genoutcome(nobs, alpha2, beta21, c(0,0), c(0,0), metx, myx2, myx3)
    probprog <- genoutcome(nobs, alpha3, c(0,0), beta2, beta3, metx, myx2, myx3)
    Zprogcov <- cbind(myx2, myx3)
  }

  # Now we construct prob with both prog and pred features
  myprob1 <- myprob2 <- matrix( 0, nrow = nobs, ncol = 3)

  for (i in 1:nobs){
    myprob1[i,] <- prob1[i,2:4]*probprog[i,2:4]/sum(prob1[i,2:4]*probprog[i,2:4])
    myprob2[i,] <- prob2[i,2:4]*probprog[i,2:4]/sum(prob2[i,2:4]*probprog[i,2:4])
  }

  myprob <- list(myprob1, myprob2)
  trtsgn <- rep(c(1,2), nobs/2)

  mytot <- array(0, dim = c(nobs, 3, nset))
  myoutot <- matrix(0, nrow = nobs, ncol = nset)
  for(i in 1:nset){
    myy <- matrix(0, nrow = nobs, ncol = 3)
    myyout <- matrix(0, nrow = nobs, ncol = 1)
    for(k in 1:nobs){
      trtemp <- trtsgn[k];
      if(trtemp == 1){
        myy[k,1:3] <- t(rmultinom(n = 1, size = 1, prob = myprob[[1]][k,]))
        }
      if(trtemp == 2){
        myy[k,1:3] <- t(rmultinom(n = 1, size = 1, prob = myprob[[2]][k,]))
        }
      myyout[k] <- match(1,myy[k,])
      trtemp <- NULL
    }
    mytot[,,i] <- myy
    myoutot[,i] <- myyout
  }

  # Aggregate biomarkers
  Xpredcov <- t(genenorm[-c(91,92),])
  cont <- 1
  if(predscen == 2){
    repeat{
      cont = cont + 1
      Xpredcov <- cbind(Xpredcov, rnorm(n = nobs, mean = 0, sd = 1))
      if(cont == nnoise){
        break
      }
    }
  }
  biom <- cbind(Xpredcov, Zprogcov)

  # RETURN
  fulldata <- list(
    Y = myy,
    Yord = myoutot,
    treatment = trtsgn,
    cov = biom,
    prob = myprob
  )
  return(fulldata)
}
