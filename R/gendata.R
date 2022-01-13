##################################################################################
## Data generating mechanism for simulation scenarios
## It follows the code by Junsheng Ma. I just re-arranged it
## A total of 92 genes are selected for 152 simuluated patients               ####
## Also simulate 100 realizations data sets                                   ####
##################################################################################

# genoutcome
#
# Generates outcome variables. The response is an ordinal outcome variable.
#
# Ma, J., Stingo, F. C., & Hobbs, B. P. (2019). Bayesian personalized
# treatment selection strategies that integrate predictive with prognostic
# determinants. \emph{Biometrical Journal}, \strong{61}(4), 902-917.
# \url{https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700323}
#
#
# @param nobs number of observations
# @param alpha intercept
# @param beta1 coefficient for the predictive covariates. More precisely it is
# a coefficient for the combination of the pca to generate the exponetial.
# @param beta2 coefficient for the first prognostic covariate (x2)
# @param beta3 coefficient for the second prognostic covariate (x3)
# @param metx combination of the pca to generate the exponetial
# @param x2 first prognostic covariate (x2)
# @param x3 second prognostic covariate (x3)
# @return a \eqn{\code{nobs}\times 4} matrix. The first column contains the generated
# outcome variable, the other three the probabilities for the categories.

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

# getdata
#
# It is just an helper function.
# Using load() or data() I had to call the dataframe by the name of the object.
# It caused the devtools::check() to return this note:
# checking R code for possible problems ... NOTE
# genmech: no visible binding for global variable ‘simupats’
# Undefined global functions or variables:
#   simupats
#
#  It is just a shortcut to get things done. I could have worked on the Lazyload
#  See https://stackoverflow.com/questions/30951204/load-dataset-from-r-package-using-data-assign-it-directly-to-a-variable
# @param ... literal character strings or names
#

getdata <- function(...)
{
  e <- new.env()
  name <- data(..., envir = e)[1]
  e[[name]]
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
#' @references
#' Ma, J., Stingo, F. C., & Hobbs, B. P. (2016). Bayesian predictive modeling
#' for genomic based personalized treatment selection. \emph{Biometrics},
#' \strong{72}(2), 575-583.
#' \url{https://onlinelibrary.wiley.com/doi/full/10.1111/biom.12448}
#'
#' Ma, J., Stingo, F. C., & Hobbs, B. P. (2019). Bayesian personalized
#' treatment selection strategies that integrate predictive with prognostic
#' determinants. \emph{Biometrical Journal}, \strong{61}(4), 902-917.
#' \url{https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700323}
#'
#' @param npred number of predictive covariates used to generate the outcome
# @param mydata input dataset file. It must be a dataset with observation on the rows and
#' variables on the columns. Covariates should not be scaled.
#' @param progscen Prognostic covariates option:
#'   0 - No prognostic biomarkers.
#'   1 - Prognostic biomarkers are considered in the original scale (default).
#'   2 - Prognostic biomarkers are transformed.
#' @param predscen Predictive covariates option:
#'   1 - \code{npred} Predictive biomarkers are considered to generate outcomes (default)
#'   2 - \code{nnoise} noisy std normals are added to the design matrix of predictive
#'   covariates in addition to 10 the predictive biomarkers considered to
#'   generate outcomes
#' @param nnoise number of noisy covariates added to predictive biomarkers
#' @param nset number of replicated scenarios generated
#' @param save logical. if TRUE the function save the results in a .rda file.
#' Default is FALSE.
#' @param filename Name given to the file is results are save in .rda file.
#' Default is ``myscenario.rda''
#' @return a list of 5 elements.
#'   \itemize{
#'   \item 1 - List of \code{nset} outcome variables.
#'   \item 2 - List of \code{nset} outcome variables in ordinal notation.
#'   \item 3 - A vector of assigned treatments.
#'   \item 4 - A matrix of all biomarkers.
#'   \item 5 - A list of outcome probabilities (one element for each treatment).
#' }
#' @export

genmech <- function(npred = 10, progscen = 1,
                    predscen = 1, nnoise = 15, nset = 30, save = FALSE,
                    filename = "myscenario"){

  mydata <- getdata("simupats")
  genenorm <- scale(as.matrix(mydata))
  mypca <- prcomp(genenorm[,c(1:npred)])

  nobs <- length(mypca$x[,1])

  # Predictive Markers
  ## use the combination of the pca to generate the exponential
  metx1 <- scale(mypca$x[,1]) + scale(mypca$x[,2]) + 0.50
  metx <- sign(metx1)*(sign(metx1)*metx1)^(0.2) * 0.45
  ## for the noisy scenario

  # Prognostic Markers
  if(progscen == 0){
    ## no prognostic covariates
    z2 <- rep(0, nobs)
    z3 <- rep(0, nobs)
  }
  if(progscen == 1){
    ## original scale
    z2 <- genenorm[,91]
    z3 <- genenorm[,92]
  }
  if(progscen == 2){
    ## transformation
    z2 <- genenorm[,91]
    z3 <- genenorm[,92]
    z2 <- sign(z2)*(sign(z2)*z2)^(0.5)
    z3 <- sign(z3)*(sign(z3)*z3)^(0.2)
  }

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

  #probabilities for treatment 1
  prob1 <- genoutcome(nobs, alpha1, beta11 ,c(0,0), c(0,0), metx, z2, z3)
  #probabilities for treatment 2
  prob2 <- genoutcome(nobs, alpha2, beta21, c(0,0), c(0,0), metx, z2, z3)
  probprog <- genoutcome(nobs, alpha3, c(0,0), beta2, beta3, metx, z2, z3)
  Zprogcov <- cbind(z2, z3)

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
  Xpredcov <- genenorm[,c(1:npred)]
  cont <- 0
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
  if(save == FALSE){
    fulldata <- list(
    Y = myy,
    Yord = myoutot,
    treatment = trtsgn,
    cov = biom,
    prob = myprob)
    return(fulldata)
  } else {
    #this is mainly done for compatibility with Ma's script. In this way
    #comparison can be run smoothly
    mydata <- Xpredcov#genenorm[,c(1:npred)]
    orgz <- cbind(genenorm[,91], genenorm[,92])
    myz2 <- z2
    myz3 <- z3
    newz <- cbind(rnorm(n = nobs, mean = 0, sd = 1), rnorm(n = nobs, mean = 0,
                                                           sd = 3))
    myfile <- paste0("data/", filename, ".rda")
    save(myoutot, mytot, mydata, trtsgn, myprob, orgz, myz2, myz3, newz,
         file = myfile)
  }
}

pred_sample <- function(p = 25, o = .80){
  ov <- floor(p*o)
  all <- c(1:p)
  comm <- sample(x = all, size = ov, replace = F)
  all <- all[-comm]
  res <- length(all)
  all <- sample(all, res)

  g1 <- all[1:floor(res/2)]
  g2 <- all[(floor(res/2)+1):res]
  g1 <- c(g1, comm)
  g2 <- c(g2, comm)
  return(list(g1 = sort(g1), g2=sort(g2)))
}

tt <- function(pred, prog){
  genenorm <- pred
  mypca <- prcomp(genenorm)

  nobs <- length(mypca$x[,1])

  # Predictive Markers
  ## use the combination of the pca to generate the exponential
  metx1 <- scale(mypca$x[,1]) + scale(mypca$x[,2]) + 0.50
  metx <- sign(metx1)*(sign(metx1)*metx1)^(0.2) * 0.45
  ## for the noisy scenario

  # Prognostic Markers
  ## transformation
  z2 <- prog[,1]
  z3 <- prog[,2]
  z2 <- sign(z2)*(sign(z2)*z2)^(0.5)
  z3 <- sign(z3)*(sign(z3)*z3)^(0.2)

  # pmts probabilities for treatment 1
  alpha1 <- c(-0.5, -1)
  beta11 <- c(2, 2.6)
  # pmts probabilities for treatment 2
  alpha2 <- c(0.7, -1)
  beta21 <- c(-1, -3)
  # pmts probabilities with prognostic only
  alpha3 <- c(1, -0.5); beta2 <- c(1, 0.5); beta3 <- c(0.7,1)

  #probabilities for treatment 1
  prob1 <- genoutcome(nobs, alpha1, beta11 ,c(0,0), c(0,0), metx, z2, z3)
  #probabilities for treatment 2
  prob2 <- genoutcome(nobs, alpha2, beta21, c(0,0), c(0,0), metx, z2, z3)
  probprog <- genoutcome(nobs, alpha3, c(0,0), beta2, beta3, metx, z2, z3)
  prog <- cbind(z2, z3)

  # Now we construct prob with both prog and pred features
  myprob1 <- myprob2 <- matrix( 0, nrow = nobs, ncol = 3)

  for (i in 1:nobs){
    myprob1[i,] <- prob1[i,2:4]*probprog[i,2:4]/sum(prob1[i,2:4]*probprog[i,2:4])
    myprob2[i,] <- prob2[i,2:4]*probprog[i,2:4]/sum(prob2[i,2:4]*probprog[i,2:4])
  }

  myprob <- list(myprob1, myprob2)
  trtsgn <- rep(c(1,2), nobs/2)

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

  return(list(myoutot = myyout, mytot = myy, pred = pred, prog = prog,
         trtsgn = trtsgn, myprob = myprob))
}

#' genmech alternative
#' @param npred number of predictive covariates used to generate the outcome
#' @param nset number of replicated scenarios generated
#' @param overlap proportion of predictors used to generate the response in
#' both the train and the validation set
#' @param data dataset to be used as input for predictive and prognostic
#'
#' @export
genmech_alt <- function(npred = 10, nset = 30, overlap = 0.8, data = "simupats"){
  set.seed(121)
  if(data == "simupats"){
    EXT = 0
  }
  mydata <- getdata(data)
  genenorm <- scale(as.matrix(mydata))

  if(EXT == 0){
    pred <- genenorm[,c(1:npred)]#restituisco questi, ma riordinati
    prog <- genenorm[,c(91:92)]#restituisco questi, ma riordinati
    if(npred > 90) stop("Using the simupats dataset the maximum number of predictive covariates is 90.")
  } else {
    pred <- genenorm[,c(1:npred)]#restituisco questi, ma riordinati
    prog <- genenorm[,c(51:53)]#restituisco questi, ma riordinati
    if(npred > 50) stop("Using the simupats_ext dataset the maximum number of predictive covariates is 50.")
  }
  groups <- pred_sample(p = npred, o = overlap)
  if(EXT == 0){
    id_train <- c(1:124)
    id_test <- c(125:152)
  } else {
    id_train <- c(1:200)
    id_test <- c(201:228)
  }

  yord <- ymat <- predmk <- progmk <- trtsgn <- prob <- vector("list", length = nset)

  for(i in 1:nset){
    train <- tt(pred = pred[id_train,groups$g1], prog[id_train,])
    test <- tt(pred[id_test,groups$g2], prog[id_test,])

    yord[[i]] <- rbind(train$myoutot, test$myoutot)
    ymat[[i]] <- rbind(train$mytot, test$mytot)
    predmk[[i]] <- pred[c(id_train, id_test),]
    progmk[[i]] <- prog[c(id_train, id_test),]
    trtsgn[[i]] <- c(train$trtsgn, test$trtsgn)#[c(id_train[,i], id_test[,i])]
    prob[[i]] <- list(rbind(train$myprob[[1]], test$myprob[[1]]),
                 rbind(train$myprob[[2]], test$myprob[[2]]))
  }

  return(list(yord = yord, ymat = ymat, pred = predmk, prog = progmk, trtsgn = trtsgn,
              prob = prob))
}
