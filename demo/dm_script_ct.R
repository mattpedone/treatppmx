rm(list=ls())

#devtools::load_all()
library(mvtnorm)
library(mcclust)
library(mclust)
library(coda)
library(mcclust.ext)
library(pROC)
library(multiROC)
require(ggplot2)
library(dplyr)
library(lessR)

Rcpp::sourceCpp(file = "src/utils_ct.cpp")
Rcpp::sourceCpp(file = "src/dm_ppmx_ct.cpp")
source(file = "R/dm_ppmx_ct.R")
source(file = "R/rppmx.R")
load(file = "data/SimuOutsce2.rda")
#mydata contiene le covariate continue
#myoutot contiene la risposta in ciascuna delle 100 repliche
# è la versione "sintetica" di mytot
#myprob contiene le probabilità di ciascun livello della risposta nei due trattamenti
#per ciascun soggetto
#newx sono due var normali di media 0 e var 1 e 3 (non usate??)
#orgx contiene le due prognostiche non trasformate
#myx2 myx3 sono le due prognostiche trasformate
#trtsgn sono è il trattamento ottimale
#set.seed(121)

#trtsgn <- sample(1:3, 152, replace = TRUE)
correct_insample <- accuracy <- c()
#for(k in 1:10){
k=1
X <- data.frame(t(mydata))[, 1:50]
#Z <- data.frame(t(mydata))[, 91:92]
Z <- data.frame(cbind(myx2, myx3))#data.frame(orgx)#
Y <- mytot[,,1]#mytot[,,sample(1:100, 1)]#
idx <- sample(1:nrow(Y), 5, replace = F)#
wk <- c(0, .4, 1)
treattest <- trtsgn[idx]#treat[idx]
trt <- trtsgn[-idx]
Xtest <- X[idx,]
X <- X[-idx,]
Ytest <- Y[idx,]
Y <- Y[-idx,]
Ztest <- Z[idx,]
Z <- Z[-idx,]
nobs <- nrow(Y)

modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y))
modelpriors$hP0_L0 <- diag(1, ncol(Y))
modelpriors$hP0_nu0 <- ncol(Y) + 2
modelpriors$hP0_V0 <- diag(1, ncol(Y))

alpha_DP <- 2
n_aux <- 5
vec_par <- c(0.0, 10.0, .5, 1.0, 2.0, 2.0, 0.1)
#double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
iterations <- 50000
burnin <- 25000
thinning <- 10

uh = F
nout <- (iterations-burnin)/thinning
time_ppm <- system.time(
  out_ppmx <- my_dm_ppmx_ct(y = Y, X = X, Xpred = Xtest,
                        z = Z, zpred = Ztest, asstreat = trt, #treatment,
                        alpha=alpha_DP, CC = n_aux, reuse = 1,
                        PPMx = 1, similarity = 2, consim=1, calibration=2,
                        similparam = vec_par, modelpriors, update_hierarchy = uh,
                        iter=iterations,burn=burnin,thin=thinning, hsp = F))

sum(apply(round(apply(out_ppmx$isypred, c(1,2), mean))==Y, 1, sum)==3)/nobs
#sum(apply(round(apply(out_ppmx$pi_out, c(1,2), mean))==Y, 1, sum)==3)/nobs

label <- out_ppmx$label
truelabel <- out_ppmx$asstreat
#num_treat <- out_ppmx$num_treat
#ntreat <- 2
#for(t in 1:ntreat){
#  currtreat <- treatment
#  mat <- label[[currtreat]]
#  mat <- mat[1:num_treat[[currtreat]],]
#  cls <- as.matrix(t(mat))
#  psm <- comp.psm(cls)
#  mc_b <- minbinder.ext(psm)
#  mc_vi <- minVI(psm)
#  ari_b[t] <- adjustedRandIndex(mc_b$cl, data$labeltrain[which(data$treat == currtreat)])
#  ari_vi[t] <- adjustedRandIndex(mc_vi$cl, data$labeltrain[which(data$treat == currtreat)])
#  ess[t] <- effectiveSize(output$nclu[t,])
#}
A0 <- out_ppmx$pipred
A0 <- array(unlist(A0), dim = c(dim(A0[[1]]), length(A0)))
A0 <- apply(A0, c(1, 2, 3), mean)
Al <- list()
for(t in 1:dim(A0)[3]){
  Al[[t]] <- A0[,,t]
}
#cat("og: ", orig_names, "\n")
names(Al) <- c("1", "2")#orig_names
#Al <- Al[order(names(Al))]
print(Al)
treatpred <- c()
if(length(idx)==1){
  for(i in 1:length(idx)){
    scoreA <- sum(Al[[1]]*wk)
    scoreB <- sum(Al[[2]]*wk)

    if(scoreA > scoreB){
      treatpred[i] <- "1"
    } else {
      treatpred[i] <- "2"
    }
  }

  print(treatpred)
  print(treattest)

  accuracy[k] <- treatpred==treattest

} else {
  for(i in 1:length(idx)){
    scoreA <- sum(Al[[1]][i,]*wk)
    scoreB <- sum(Al[[2]][i,]*wk)
    #print(scoreA);print(scoreB);
    if(scoreA > scoreB){
      treatpred[i] <- "1"
    } else {
      treatpred[i] <- "2"
    }
  }

  print(treatpred)
  print(treattest)

  accuracy[k] <- sum(diag(table(treatpred, treattest)))/length(idx)

}
#}

#correct_insample
accuracy

#A <- array(unlist(out_ppmx$eta), dim = c(dim(out_ppmx$eta[[1]]), nout, ntreat))
#idx <- which(ncol(A[,,,2])==9)
#apply(A[,,,1], c(1, 2), mean)[,1:10]
#apply(A[,,,2], c(1, 2), mean)[,1:10]
