#rm(list=ls())
#
#library(Rcpp)
#library(RcppArmadillo)
#library(mvtnorm)
#library(mcclust)
#library(mclust)
#library(coda)
#library(mcclust.ext)
##require(ggplot2)
##library(reshape2)
#
#Rcpp::sourceCpp(file = "src/utils_ct.cpp")
#Rcpp::sourceCpp(file = "src/dm_ppmx_ct.cpp")
#source(file = "R/rppmx.R")
#source(file = "R/dm_ppmx_ct.R")
#
#load("data/SimuOutsce2.rda")
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

devtools::load_all()
K <- 10 #repliche
risultati <- array(0, dim = c(12, 9, K))

vecadp <- c(1, 2, 10)
vecagow <- c(1, 5, 10, 20)
idxsc <- 1


for(alphagow in 1:4){
  for(alphadp in 1:3){
    for(k in 1:K){
      X <- data.frame(t(mydata))[, -c(11:92)]#data.frame(mydata)#
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
      alpha_gow <- vecagow[alphagow]
      n_aux <- 5
      vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
      #double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
      iterations <- 75000; burnin <- 25000; thinning <- 10

      nout <- (iterations-burnin)/thinning
      time_ppmx <- system.time(
        out_ppmx <- my_dm_ppmx_ct(y = Y, X = X, Xpred = Xtest,
                                  z = Z, zpred = Ztest, asstreat = trt, #treatment,
                                  alpha = alpha_DP, CC = n_aux, reuse = 1,
                                  PPMx = 1, similarity = 3, consim = 1,  gowtot = 1,
                                  alphagow = alpha_gow, calibration = 2, coardegree = 1,
                                  similparam = vec_par, modelpriors, update_hierarchy = T,
                                  iter = iterations, burn = burnin, thin = thinning, hsp = T))
      risultati[idxsc,9,k] <- as.double(time_ppmx[3])
      # Posterior clustering ----
      num_treat <- table(trt)
      risultati[idxsc,1,k] <- mean(apply(out_ppmx$label[[1]], 2, max))
      cls <- t(as.matrix(out_ppmx$label[[1]]))[,c(1:num_treat[1])]
      psm <- comp.psm(cls)
      mc_vi <- minVI(psm); risultati[idxsc,3,k] <- max(mc_vi$cl)
      risultati[idxsc,2,k] <- mean(apply(out_ppmx$label[[2]], 2, max))
      cls2 <- t(as.matrix(out_ppmx$label[[2]]))[,c(1:num_treat[2])]
      psm2 <- comp.psm(cls2)
      mc_b2 <- minbinder.ext(psm2); max(mc_b2$cl)
      mc_vi2 <- minVI(psm2); risultati[idxsc,4,k] <- max(mc_vi2$cl)
      # In sample prediction (goodness-of-fit) ----
      # overall
      risultati[idxsc,7,k] <- sum(apply(round(apply(out_ppmx$isypred, c(1,2), mean))==Y, 1, sum)==3)/nobs
      # by treatment
      risultati[idxsc,5,k] <- sum(apply(round(apply(out_ppmx$isypred[which(trt == 1),,], c(1,2), mean))==Y[which(trt == 1),], 1, sum)==3)/sum((trt == 1))
      risultati[idxsc,6,k] <- sum(apply(round(apply(out_ppmx$isypred[which(trt == 2),,], c(1,2), mean))==Y[which(trt == 2),], 1, sum)==3)/sum((trt == 2))
      #posterior predictive probabilities ----
      A0 <- apply(out_ppmx$ypred, c(1,2,3), mean);#A0
      #treatmente prediction with utility function ----
      optrt <- as.numeric(myprob[[2]][idx,]%*%wk > myprob[[1]][idx,]%*%wk)+1
      predtrt <- as.numeric(A0[,,2]%*%wk > A0[,,1]%*%wk)+1
      risultati[idxsc,8,k] <- sum(optrt==predtrt)/length(predtrt)
    }
    idxsc <- idxsc+1
  }
}

#calcola indici paper
#fai script x confronto

apply(risultati, c(1,2), mean)
sqrt(apply(risultati, c(1,2), var))

#[,1] [,2]      [,3]      [,4] [,5]      [,6]    [,7]
#[1,]    1    1 0.8000000 0.7400000 0.77 0.3461538 143.718
#[2,]    1    1 0.8431373 0.8163265 0.83 0.3269231 156.441
#[3,]    1    1 1.0000000 1.0000000 1.00 0.6346154 235.398
#[4,]    1    1 0.7800000 0.7600000 0.77 0.3076923 148.672
#[5,]    1    1 0.7843137 0.7142857 0.75 0.2692308 157.999
#[6,]    1    1 1.0000000 0.9795918 0.99 0.3653846 226.185
#[7,]    1    1 0.6875000 0.7307692 0.71 0.3461538 140.929
#[8,]    1    1 0.9000000 0.7000000 0.80 0.3846154 153.299
#[9,]    1    1 1.0000000 0.9583333 0.98 0.5000000 230.045
#[10,]    1    1 0.7800000 0.7200000 0.75 0.3653846 150.401
#[11,]    1    1 0.8723404 0.7547170 0.81 0.4038462 160.246
#[12,]    1    1 0.9795918 1.0000000 0.99 0.7115385 240.185
#[13,]    1    1 0.7021277 0.7547170 0.73 0.2884615 147.772
#[14,]    1    1 0.8163265 0.7058824 0.76 0.2884615 157.651
#[15,]    1    1 1.0000000 0.9607843 0.98 0.3076923 229.274
#[16,]    1    1 0.9387755 0.8627451 0.90 0.3653846 165.269
#[17,]    1    1 0.9411765 0.9591837 0.95 0.4423077 174.388
#[18,]    1    1 1.0000000 1.0000000 1.00 0.4807692 256.294
#[19,]    1    1 0.9038462 0.8333333 0.87 0.2884615 170.946
#[20,]    1    1 0.9803922 0.9387755 0.96 0.4615385 182.306
#[21,]    1    1 1.0000000 1.0000000 1.00 0.4807692 279.404
#[22,]    1    1 0.9019608 0.9795918 0.94 0.5769231 195.712
#[23,]    1    1 0.9782609 0.9814815 0.98 0.4807692 210.554
#[24,]    1    1 1.0000000 1.0000000 1.00 0.6346154 305.669
#[25,]    1    1 1.0000000 0.9565217 0.98 0.7500000 232.961
#[26,]    1    1 1.0000000 1.0000000 1.00 0.6346154 256.021
#[27,]    1    1 1.0000000 1.0000000 1.00 0.6346154 328.437
#[28,]    1    1 1.0000000 1.0000000 1.00 0.6538462 275.171
#[29,]    1    1 1.0000000 1.0000000 1.00 0.5576923 305.237
#[30,]    2    2 1.0000000 1.0000000 1.00 0.4615385 356.561#
