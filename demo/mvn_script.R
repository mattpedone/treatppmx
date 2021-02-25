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
library(coda)

sourceCpp("src/mvn_ppmx.cpp")
source("R/mvn_ppmx.R")
source("R/rppmx.R")
#set.seed(24)
K=2
myppmx <- gcd(n=100, concov = 10, K, alpha = 5)
Y <- myppmx$y
X <- myppmx$X

colors <- c("#ebb678", "#1979a9", "#e07b39", "#69bdd2", "#80391e", "#cce7e8", 
            "#1c100b", "#042f66", "#44bcd8")
colors <- colors[myppmx$label]
if(K==2){plot(Y, pch = 16, col = colors)}
if(K==3){scatterplot3d(Y, pch = 16, color = colors, grid=TRUE, box=FALSE)}

modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y))
modelpriors$hP0_L0 <- diag(10, ncol(Y))
modelpriors$hP0_nu0 <- nrow(Y) + 2
modelpriors$hP0_V0 <- diag(10, ncol(Y))

mytime <- system.time(
  out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 3, reuse = 1, PPMx = 1, 
                     similarity = 1, consim=1, calibration=2, 
                     similparam=c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                     modelpriors, mhtune=c(0.5, 0.5), iter=10000, burn=5000, 
                     thin=10))

postquant()
myppmx$nclus

apply(out$mu, c(1, 2), mean)
out$mu[,,250]
myppmx$possmean

colors <- c("#ebb678", "#1979a9", "#e07b39", "#69bdd2", "#80391e", "#cce7e8", 
            "#1c100b", "#042f66", "#44bcd8")
colors <- colors[postquant(lab = T)$lab]
if(K==2){plot(Y, pch = 16, col = colors)}
if(K==3){scatterplot3d(Y, pch = 16, color = colors, grid=TRUE, box=FALSE)}


