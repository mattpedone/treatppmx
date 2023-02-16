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
    temp <- beta1*(metx[i])^3 + beta2*x2[i] + beta3*x3[i]
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
#' Generates the data for the simulations scenarios 1a and 1b reported in the paper.
#' It follows the strategy designed by Ma et al. (2019).
#' In particular, the outcome is a categorical variable representing \eqn{K=3} benefit-increasing levels.
#' Patients (\eqn{n=152}) are assigned to \eqn{T=2} competing treatments.
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
#' @param npred number of \eqn{Q} predictive covariates used to generate the outcome
#' @param progscen Prognostic covariates option:
#'   0 - No prognostic biomarkers;
#'   1 - Prognostic biomarkers are considered in the original scale (default);
#'   2 - Prognostic biomarkers are transformed
#' @param predscen Predictive covariates option:
#'   1 - \code{npred} Predictive biomarkers are considered to generate outcomes (default);
#'   2 - \code{nnoise} noisy std normals are added to the design matrix of predictive
#'   covariates in addition to \code{npred} predictive biomarkers considered to
#'   generate outcomes
#' @param nnoise number of noisy covariates added to predictive biomarkers
#' @param nset number of replicated scenarios generated
#' @param save logical. if TRUE the function save the results in a \code{.rda} file.
#' Default is FALSE
#' @param filename Name given to the file is results are saved in \code{.rda} file
#' Default is \code{myscenario.rda}
#' @return a list of 5 elements.
#'   \itemize{
#'   \item \code{Y}: \eqn{n\times K \times}\code{nset} array. It contains \code{nset} matrices that store the outcome in the form of a multinomial experiment
#'   \item \code{Yord}: \eqn{n\times}\code{nset} matrix. It contains the ordinal outcome for each replicated dataset
#'   \item \code{treatment}: \eqn{n-}dimensional vector. It contains the treatment assigned
#'   \item \code{cov}: \eqn{n\times (Q+2)} matrix. It contains all the biomarkers. The last two columns are the two predictive biomarkers
#'   \item \code{prob}: List of \eqn{T}. Each element is a \eqn{n\times K} matrix containing the response probabilities
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
    Y = mytot,
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
  beta11 <- c(1.5, 2)
  # pmts probabilities for treatment 2
  alpha2 <- c(0.7, -1)
  beta21 <- c(-.5, -1.0)
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


#' genmech_het
#'
#' Generates the data for the simulations scenarios reported in the paper.
#' It follows the strategy designed by Ma et al. (2019).
#' In particular, the outcome is a categorical variable representing \eqn{K=3} benefit-increasing levels.
#' Patients (\eqn{n=152}) are assigned to \eqn{T=2} competing treatments.
#' @param npred number of \eqn{Q} predictive covariates used to generate the outcome
#' @param nset number of replicated scenarios generated
#' @param overlap proportion of predictors used to generate the response in
#' both the train and the validation set
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
#' @return a list of 6 elements.
#'   \itemize{
#'   \item \code{yord}: List of \code{nset}. Each element is a \code{n-}dimensional vector of the ordinal outcome
#'   \item \code{ymat}: List of \code{nset}. Each element is a \eqn{n\times K} matrix containing the ordinal outcome in the form of a Multinomial experiment
#'   \item \code{pred}: List of \code{nset}. Each element is a \eqn{n\times Q} matrix containing the predictive biomarkers
#'   \item \code{pred}: List of \code{nset}. Each element is a \eqn{n\times 2} matrix containing the prognostic biomarkers
#'   \item \code{trtsgn}: List of \code{nset}. Each element is a \code{n-}dimensional vector of the treatment assigned
#'   \item \code{prob}: List of \code{nset}. Each element of the list is a list of \eqn{T} \eqn{n\times K} matrices containing the response probabilities
#' }
#' @export

genmech_het <- function(npred = 10, nset = 30, overlap = 0.8){
  set.seed(121)
  mydata <- getdata("simupats")

  #for(i in 1:nset){
  genenorm <- scale(as.matrix(mydata))
  #genenorm <- genenorm[sample(nrow(genenorm)),]
  if(npred > 90) stop("Using the simupats dataset the maximum number of predictive covariates is 90.")
  pred <- genenorm[,c(1:npred)]#restituisco questi, ma riordinati
  prog <- genenorm[,c(91:92)]#restituisco questi, ma riordinati

  groups <- pred_sample(p = npred, o = overlap)
  id_train <- c(1:124)
  id_test <- c(125:152)

  yord <- ymat <- predmk <- progmk <- trtsgn <- prob <- vector("list", length = nset)

  for(i in 1:nset){
    train <- tt(pred = pred[id_train,groups$g1], prog = prog[id_train,])
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

#' genmech clustering
#'
#' Generates the data for the simulations scenarios S1 reported in the Supplementary Material.
#' It follows the strategy designed by Ma et al. (2019).
#' In particular, the outcome is a categorical variable representing \eqn{K=3} benefit-increasing levels.
#' Patients (\eqn{n}) are assigned to \eqn{T=2} competing treatments.
#' Moreover, the clustering is known and fixed. We generate predictive biomarkers as in Argiento et al. (2022).
#'
#' @param npred number of \eqn{Q} predictive covariates used to generate the outcome
#' @param n number of observations
#' @param nnoise number of noisy variables
# @param nset number of replicated scenarios generated
#'
#' @references
#' Argiento, R., Corradin, R., and Guglielmi, A. (2022). A Bayesian nonparametric model
#' for covariate driven clustering: improved insights of blood donors data. \emph{Unpublished
#' Manuscript}.
#'
#' Ma, J., Stingo, F. C., & Hobbs, B. P. (2019). Bayesian personalized
#' treatment selection strategies that integrate predictive with prognostic
#' determinants. \emph{Biometrical Journal}, \strong{61}(4), 902-917.
#' \url{https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700323}
#'
#' @return a list of 8 elements.
#'   \itemize{
#'   \item \code{Y}: \eqn{n\times K \times}\code{nset} array. It contains \code{nset} matrices that store the outcome in the form of a multinomial experiment
#'   \item \code{Yord}: \eqn{n\times}\code{nset} matrix. It contains the ordinal outcome for each replicated dataset
#'   \item \code{treatment}: \eqn{n-}dimensional vector. It contains the treatment assigned
#'   \item \code{pred}: \eqn{n\times}\eqn{(}\code{npred}\eqn{+}\code{nnoise}\eqn{)} matrix. It contains the predictive biomarkers (only the first \code{npred} are effectively used to generate the response)
#'   \item \code{pred}: \eqn{n\times 2} matrix. It contains the prognostic biomarkers
#'   \item \code{clu1}: contains the cluster label for patients assigned to Treatment 1
#'   \item \code{clu2}: contains the cluster label for patients assigned to Treatment 2
#'   \item \code{prob}: List of \eqn{T}. Each element is a \eqn{n\times K} matrix containing the response probabilities
#' }
#' @export

genmech_clu <- function(npred = 4, n = 200, nnoise = 7){#, nset = 50){

  N <- nobs <- n

  ##Variable to store the samples from the mixture distribution
  mydata <- matrix(NA, N, npred)

  U <- c(rep(1, 75), rep(2, 75), rep(3, 50))
  U <- U[sample(1:N, N, replace = TRUE)]
  for(i in 1:N){
    if(U[i]==1){
      mydata[i, 1] <- rnorm(1, -3, .5)
      mydata[i, 2] <- rnorm(1, 3, .5)
      mydata[i, 3] <- rbinom(1, 1, .1)
      mydata[i, 4] <- rbinom(1, 1, .1)
    }
    if(U[i]==2){
      mydata[i, 1] <- rnorm(1, 0, .5)
      mydata[i, 2] <- rnorm(1, 0, .5)
      mydata[i, 3] <- rbinom(1, 1, .5)
      mydata[i, 4] <- rbinom(1, 1, .5)
    }
    if(U[i]==3){
      mydata[i, 1] <- rnorm(1, 3, .5)
      mydata[i, 2] <- rnorm(1, 3, .5)
      mydata[i, 3] <- rbinom(1, 1, .9)
      mydata[i, 4] <- rbinom(1, 1, .9)
    }
  }

  genenorm <- cbind(scale(as.matrix(mydata[,1:2])), mydata[,3:4])
  #genenorm <- mydata
  mypca <- prcomp(genenorm[,c(1:npred)])

  # Predictive Markers
  ## use the combination of the pca to generate the exponential
  metx1 <- scale(mypca$x[,1]) + scale(mypca$x[,2]) + 0.50
  metx <- sign(metx1)*(sign(metx1)*metx1)^(0.2) * 0.45
  ## for the noisy scenario

  # Prognostic Markers
  z2 <- rnorm(nobs)
  z3 <- rnorm(nobs)

  z2 <- sign(z2)*(sign(z2)*z2)^(0.5)
  z3 <- sign(z3)*(sign(z3)*z3)^(0.2)

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

  #mytot <- matrix(0, nobs, 3)
  #myoutot <- rep(0, nobs)
  #for(i in 1:nset){
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
    #mytot[k,] <- myy
    #myoutot[,i] <- myyout
  #}

  # Aggregate biomarkers
  Xpredcov <- genenorm[,c(1:npred)]
  if(nnoise > 0){
    cont <- 0
    repeat{
      cont = cont + 1
      Xpredcov <- cbind(Xpredcov, rnorm(n = nobs, mean = 0, sd = 1))
      if(cont == nnoise){
        break
      }
    }
  }

  # RETURN
  fulldata <- list(
      Y = myy, #mytot,
      Yord = myyout, #myoutot,
      treatment = trtsgn,
      pred = Xpredcov,
      prog = Zprogcov,
      clu1 = U[which((1:length(U))%%2==1)],
      clu2 = U[which((1:length(U))%%2==0)],
      prob = myprob)
    return(fulldata)
}

#' genmech clustering 2
#'
#' Generates the data for the simulations scenarios S3 reported in the Supplementary Material.
#' It follows the strategy designed by Ma et al. (2019).
#' In particular, the outcome is a categorical variable representing \eqn{K=3} benefit-increasing levels.
#' Patients (\eqn{n}) are assigned to \eqn{T=2} competing treatments.
#' Moreover, the clustering is known and fixed. We generate predictive biomarkers as Scenarion 4 in Page and Quintana (2018)
#' @references
#' Ma, J., Stingo, F. C., & Hobbs, B. P. (2019). Bayesian personalized
#' treatment selection strategies that integrate predictive with prognostic
#' determinants. \emph{Biometrical Journal}, \strong{61}(4), 902-917.
#' \url{https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700323}
#'
#' Page, G. L. and Quintana, F. A. (2018). Calibrating covariate informed
#' product partition models. \emph{Statistics and Computing}, \strong{28}(5), 1009–1031.
#'
#' @param npred number of predictive covariates used to generate the outcome
#'
#' @return a list of 8 elements.
#'   \itemize{
#'   \item \code{Y}: \eqn{n\times K \times}\code{nset} array. It contains \code{nset} matrices that store the outcome in the form of a multinomial experiment
#'   \item \code{Yord}: \eqn{n\times}\code{nset} matrix. It contains the ordinal outcome for each replicated dataset
#'   \item \code{treatment}: \eqn{n-}dimensional vector. It contains the treatment assigned
#'   \item \code{pred}: \eqn{n\times}\eqn{(}\code{npred}\eqn{+}\code{nnoise}\eqn{)} matrix. It contains the predictive biomarkers (only the first \code{npred} are effectively used to generate the response)
#'   \item \code{pred}: \eqn{n\times 2} matrix. It contains the prognostic biomarkers
#'   \item \code{clu1}: contains the cluster label for patients assigned to Treatment 1
#'   \item \code{clu2}: contains the cluster label for patients assigned to Treatment 2
#'   \item \code{prob}: List of \eqn{T}. Each element is a \eqn{n\times K} matrix containing the response probabilities
#' }
#' @export

genmech_clu2 <- function(npred = 10){

  N <- nobs <- n <- 200

  ##Variable to store the samples from the mixture distribution
  mydata <- matrix(NA, N, npred)

  nc <- npred/5

  M <- matrix(0, N, nc)

  No <- matrix(rnorm(N*nc), N, nc)
  UN <- matrix(runif(N*nc, 0, 10), N, nc)
  TS <- matrix(rt(N*nc, df = 4), N, nc)
  SN <- matrix(sn::rsn(n=N*nc, xi=10, omega=1, alpha=10, tau=0,  dp=NULL), N, nc)
  pp <- runif(N)
  for(i in 1:N){
    if(pp[i]<.4){
      M[i, 1] <- rnorm(1)
      M[i, 2] <- rnorm(1)
    } else {
      M[i, 1] <- rnorm(1, 10, 2)
      M[i, 2] <- rnorm(1, 10, 2)
    }
  }

  U <- c(rep(1, 50), rep(2, 50), rep(3, 50), rep(4, 50))
  #U <- U[sample(1:N, N, replace = TRUE)]
  for(i in 1:N){
    if(U[i]==1){
      No[i, 1:nc] <- rnorm(nc, 3, 1)
      UN[i, 1:nc] <- runif(nc, -5, 0)
    }
    if(U[i]==2){
      TS[i, 1:nc] <- rt(1, df = 4) - 5
      SN[i, 1:nc] <- sn::rsn(n=nc, xi=10, omega=10, alpha=10, tau=0,  dp=NULL)
    }
    if(U[i]==3){
      for(i in 1:N){
        if(pp[i]>.4){
          M[i, 1:nc] <- rnorm(nc, -10, 2)
          #M[i, 2:nc] <- rnorm(nc, -10, 2)
        }
      }
    }
  }

  mydata <- scale(cbind(No, UN, TS, SN, M))

  #genenorm <- cbind(scale(as.matrix(mydata[,1:2])), mydata[,3:4])
  genenorm <- mydata
  mypca <- prcomp(genenorm[,c(1:npred)])

  # Predictive Markers
  ## use the combination of the pca to generate the exponential
  metx1 <- scale(mypca$x[,1]) + scale(mypca$x[,2]) + 0.50
  metx <- sign(metx1)*(sign(metx1)*metx1)^(0.2) * 0.45
  ## for the noisy scenario

  # Prognostic Markers
  z2 <- rnorm(nobs)
  z3 <- rnorm(nobs)

  z2 <- sign(z2)*(sign(z2)*z2)^(0.5)
  z3 <- sign(z3)*(sign(z3)*z3)^(0.2)

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

  #mytot <- array(0, dim = c(nobs, 3, nset))
  #myoutot <- matrix(0, nrow = nobs, ncol = nset)
  #for(i in 1:nset){
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
  #  mytot[,,i] <- myy
  #  myoutot[,i] <- myyout
  #}

  # Aggregate biomarkers
  Xpredcov <- genenorm[,c(1:npred)]

  # RETURN
  fulldata <- list(
    Y = myy, #mytot,
    Yord = myyout, #myoutot,
    treatment = trtsgn,
    pred = Xpredcov,
    prog = Zprogcov,
    clu1 = U[which((1:length(U))%%2==1)],
    clu2 = U[which((1:length(U))%%2==0)],
    prob = myprob)
  return(fulldata)
}

#' genmech clustering 3
#'
#' Generates the data for the simulations scenarios S2 reported in the Supplementary Material.
#' It follows the strategy designed by Ma et al. (2019).
#' In particular, the outcome is a categorical variable representing \eqn{K=3} benefit-increasing levels.
#' \eqn{n} patients are assigned to \eqn{T=2} competing treatments.
#' Moreover, the clustering is known and fixed. We generate predictive biomarkers as Scenarion 4 in Page and Quintana (2018)
#'
#' @references
#' Ma, J., Stingo, F. C., & Hobbs, B. P. (2019). Bayesian personalized
#' treatment selection strategies that integrate predictive with prognostic
#' determinants. \emph{Biometrical Journal}, \strong{61}(4), 902-917.
#' \url{https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700323}
#'
#' Page, G. L. and Quintana, F. A. (2018). Calibrating covariate informed
#' product partition models. \emph{Statistics and Computing}, \strong{28}(5), 1009–1031.
#'
#' @param nset number of replicated scenarios generated
#'
#' @return a list of 8 elements.
#'   \itemize{
#'   \item \code{Y}: \eqn{n\times K \times}\code{nset} array. It contains \code{nset} matrices that store the outcome in the form of a multinomial experiment
#'   \item \code{Yord}: \eqn{n\times}\code{nset} matrix. It contains the ordinal outcome for each replicated dataset
#'   \item \code{treatment}: \eqn{n-}dimensional vector. It contains the treatment assigned
#'   \item \code{pred}: \eqn{n\times}\eqn{(}\code{npred}\eqn{+}\code{nnoise}\eqn{)} matrix. It contains the predictive biomarkers (only the first \code{npred} are effectively used to generate the response)
#'   \item \code{pred}: \eqn{n\times 2} matrix. It contains the prognostic biomarkers
#'   \item \code{clu1}: contains the cluster label for patients assigned to Treatment 1
#'   \item \code{clu2}: contains the cluster label for patients assigned to Treatment 2
#'   \item \code{prob}: List of \eqn{T}. Each element is a \eqn{n\times K} matrix containing the response probabilities
#' }
#'
# @export

genmech_clu3 <- function(nset = 50){

  N <- nobs <- 200

  ##Variable to store the samples from the mixture distribution
  #x1 <- scale(rnorm(N))
  #x2 <- scale(runif(N, 0, 10))
  x1 <- rnorm(N)
  x2 <- runif(N, 0, 10)
  x3 <- rbinom(N, 1, .5)
  x4 <- sample(factor(c("a", "b", "c")), N, replace =T)

  o <- lm(x2~x1+x3+x4+x1:x3+x1:x4+x3:x4+x1:x3:x4, x=TRUE)
  mydata <- cbind(o$x, x2)

  groups <- mydata[, c(3:5)]
  U <- vector()
  for(i in 1:N){
    if(sum(groups[i,] == c(0, 0, 0)) == 3){
      U[i] <- 1
    }
    if(sum(groups[i,] == c(1, 0, 0)) == 3){
      U[i] <- 2
    }
    if(sum(groups[i,] == c(0, 1, 0)) == 3){
      U[i] <- 3
    }
    if(sum(groups[i,] == c(1, 1, 0)) == 3){
      U[i] <- 4
    }
    if(sum(groups[i,] == c(0, 0, 1)) == 3){
      U[i] <- 5
    }
    if(sum(groups[i,] == c(1, 0, 1)) == 3){
      U[i] <- 6
    }
  }

  npred <- dim(mydata)[2]
  #mydata[,c(2, 6:8, 10:12)] <- scale(mydata[,c(2, 6:8, 10:12)])
  #mydata

  #genenorm <- cbind(scale(as.matrix(mydata[,1:2])), mydata[,3:4])
  genenorm <- mydata
  #mypca <- prcomp(genenorm[,c(1:npred)])
  #
  ## Predictive Markers
  ### use the combination of the pca to generate the exponential
  #metx1 <- scale(mypca$x[,1]) + scale(mypca$x[,2]) + 0.50
  coeff <- c(-1, -2, 6, -4, 2, 5/2, 3/2, 1, -11, 3, -4, 1/2, 0.0)
  metx1 <- genenorm%*%as.vector(coeff)#+rnorm(nobs)
  metx <- sign(metx1)*(sign(metx1)*metx1)^(0.2) * 0.45
  ## for the noisy scenario

  # Prognostic Markers
  z2 <- rnorm(nobs)
  z3 <- rnorm(nobs)

  z2 <- sign(z2)*(sign(z2)*z2)^(0.5)
  z3 <- sign(z3)*(sign(z3)*z3)^(0.2)

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
  Xpredcov <- cbind(x1, x2, x3, x4)

  # RETURN
  fulldata <- list(
    Y = mytot,
    Yord = myoutot,
    treatment = trtsgn,
    pred = Xpredcov,
    prog = Zprogcov,
    clu1 = U[which((1:length(U))%%2==1)],
    clu2 = U[which((1:length(U))%%2==0)],
    prob = myprob)
  return(fulldata)
}
