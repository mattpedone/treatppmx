rm(list=ls())

library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)
library(mcclust)
library(mclust)
library(coda)
library(mcclust.ext)
#require(ggplot2)
#library(reshape2)

Rcpp::sourceCpp(file = "src/utils_ct.cpp")
Rcpp::sourceCpp(file = "src/dm_ppmx_ct.cpp")
source(file = "R/rppmx.R")
source(file = "R/dm_ppmx_ct.R")

load("data/SimuOutsce2.rda")
#data from Ma et al. Biom. J. (2017)
#mydata contiene le covariate continue
#myoutot contiene la risposta in ciascuna delle 100 repliche
# è la versione "sintetica" di mytot
#myprob contiene le probabilità di ciascun livello della risposta nei due trattamenti
#per ciascun soggetto
#newx sono due var normali di media 0 e var 1 e 3 (non usate??)
#orgx contiene le due prognostiche non trasformate
#myx2 myx3 sono le due prognostiche trasformate
#trtsgn sono i trattamenti assegnati casualmente (trial clinico)

K <- 1 #repliche
risultati <- array(0, dim = c(6, 7, K))

vecadp <- c(1, 2, 10)
idxsc <- 1

for(consim in 1:2){
  for(alphadp in 1:3){
    for(k in 1:K){
      X <- data.frame(t(mydata))[, -c(41:92)]#data.frame(mydata)#
      Z <- data.frame(cbind(myx2, myx3))#data.frame(orgx)#
      Y <- mytot[,,k]
      idx <- sort(sample(1:nrow(Y), 52, replace = F))
      wk <- c(0, 40, 100)

      treattest <- trtsgn[idx]; trt <- trtsgn[-idx]
      Xtest <- data.frame(X[idx,]); X <- data.frame(X[-idx,])
      Ytest <- data.matrix(Y[idx,]); Y <- data.matrix(Y[-idx,])
      Ztest <- data.frame(Z[idx,]); Z <- data.frame(Z[-idx,])
      nobs <- nrow(Y)

      optrt <- as.numeric(myprob[[2]][idx,]%*%wk > myprob[[1]][idx,]%*%wk)+1; #optrt

      modelpriors <- list()
      modelpriors$hP0_m0 <- rep(0, ncol(Y)); modelpriors$hP0_L0 <- diag(10, ncol(Y))
      modelpriors$hP0_nu0 <- ncol(Y) + 2; modelpriors$hP0_V0 <- diag(10, ncol(Y))

      alpha_DP <- vecadp[alphadp]
      n_aux <- 5
      vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
      #double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
      iterations <- 75000; burnin <- 25000; thinning <- 10

      nout <- (iterations-burnin)/thinning
      time_ppmx <- system.time(
        out_ppmx <- my_dm_ppmx_ct(y = Y, X = X, Xpred = Xtest,
                                  z = Z, zpred = Ztest, asstreat = trt, #treatment,
                                  alpha = alpha_DP, CC = n_aux, reuse = 1,
                                  PPMx = 1, similarity = 2, consim = consim,  gowtot = 1,
                                  alphagow = 5, calibration = 1, coardegree = 1,
                                  similparam = vec_par, modelpriors, update_hierarchy = T,
                                  iter = iterations, burn = burnin, thin = thinning, hsp = T))
      risultati[idxsc,7,k] <- as.double(time_ppmx[3])

      # Posterior clustering ----
      num_treat <- table(trt)
      cls <- t(as.matrix(out_ppmx$label[[1]]))[,c(1:num_treat[1])]
      psm <- comp.psm(cls)
      mc_vi <- minVI(psm); risultati[idxsc,1,k] <- max(mc_vi$cl)


      cls2 <- t(as.matrix(out_ppmx$label[[2]]))[,c(1:num_treat[2])]
      psm2 <- comp.psm(cls2)
      mc_b2 <- minbinder.ext(psm2); max(mc_b2$cl)
      mc_vi2 <- minVI(psm2); risultati[idxsc,2,k] <- max(mc_vi2$cl)

      # In sample prediction (goodness-of-fit) ----
      # overall
      risultati[idxsc,5,k] <- sum(apply(round(apply(out_ppmx$isypred, c(1,2), mean))==Y, 1, sum)==3)/nobs
      # by treatment
      risultati[idxsc,3,k] <- sum(apply(round(apply(out_ppmx$isypred[which(trt == 1),,], c(1,2), mean))==Y[which(trt == 1),], 1, sum)==3)/sum((trt == 1))
      risultati[idxsc,4,k] <- sum(apply(round(apply(out_ppmx$isypred[which(trt == 2),,], c(1,2), mean))==Y[which(trt == 2),], 1, sum)==3)/sum((trt == 2))

      #posterior predictive probabilities ----
      A0 <- apply(out_ppmx$ypred, c(1,2,3), mean);#A0

      #treatmente prediction with utility function ----
      optrt <- as.numeric(myprob[[2]][idx,]%*%wk > myprob[[1]][idx,]%*%wk)+1
      predtrt <- as.numeric(A0[,,2]%*%wk > A0[,,1]%*%wk)+1

      risultati[idxsc,6,k] <- sum(optrt==predtrt)/length(predtrt)
    }
    idxsc <- idxsc+1
  }

}
#calcola indici paper
#fai script x confronto

apply(risultati, c(1,2), mean)
sqrt(apply(risultati, c(1,2), var))

#[,1] [,2]      [,3]      [,4] [,5]      [,6]    [,7]
#[1,]    3    3 0.5714286 0.6666667 0.62 0.1730769 285.742
#[2,]    2    3 0.8125000 0.7500000 0.78 0.2692308 260.829
#[3,]    3    3 0.7000000 0.7200000 0.71 0.1153846 285.730
#[4,]    2    6 0.7272727 0.6250000 0.67 0.2884615 640.577
#[5,]    6    3 0.5625000 0.7692308 0.67 0.4807692 826.900
#[6,]    3    7 0.6800000 0.6800000 0.68 0.2500000 740.978#
