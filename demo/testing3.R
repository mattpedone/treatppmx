rm(list=ls())
if(!is.null(dev.list())) dev.off()
# Clear console
#cat("\014") 
#devtools::load_all()
library(Rcpp)
library(rrr)
library(dplyr)
library(tidyselect)
library(mcclust)
library(mcclust.ext)
library(mclust)

sourceCpp("src/mvn_ppmx.cpp")
source("R/mvn_ppmx.R")
source("R/gendata.R")
source("R/rppmx.R")
#set.seed(24)
n=10
d <- genera_dati(n=n)

X <- d$XX

X$X1 <- as.factor(X$X1)
X$X2 <- as.factor(X$X2)

myppmx <- ran_ppmx(X=X, similarity = 1, simparm = 1, alpha=1, m0=0, s20=1,v=2, k0=10, v0=1)
myppmx

par(mfrow=c(1, 2))
possmean <- matrix(0, myppmx$nclus, 2)
for(i in 1:myppmx$nclus){
  possmean[i,] <- mvtnorm::rmvnorm(1, c(0, 0), sigma = diag(25, 2), method="svd")
}
plot(possmean)

Y <- matrix(0, n, 2)
for(i in 1:n){
  Y[i, ] <- mvtnorm::rmvnorm(1, possmean[myppmx$label[i],], sigma = diag(.25, 2), method="svd")
}
plot(Y)

modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y))
modelpriors$hP0_L0 <- diag(10, ncol(Y))
modelpriors$hP0_nu0 <- nrow(Y) + 2
modelpriors$hP0_V0 <- diag(10, ncol(Y))

mbm <- microbenchmark::microbenchmark("nocal" =  my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 3, reuse = 1, PPMx = 1, similarity = 1, consim=1, calibration=0, 
                                                           similparam=c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                                                           modelpriors, mhtune=c(0.5, 0.5), 
                                                           iter=100,burn=0,thin=1), 
                                      "cal" =  my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 3, reuse = 1, PPMx = 1, similarity = 1, consim=1, calibration=2, 
                                                            similparam=c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                                                            modelpriors, mhtune=c(0.5, 0.5), 
                                                            iter=100,burn=0,thin=1), times = 1000)

mbm
ggplot2::autoplot(mbm)
