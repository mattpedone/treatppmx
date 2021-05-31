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
load(file = "data/SimOutScen1.RData")
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
for(k in 1:3){
X <- data.frame(t(mydata))#[, 1:4]
Z <- data.frame(orgx)#data.frame(cbind(myx2, myx3))
Y <- mytot[,,k]
idx <- sample(1:nrow(Y), 5, replace = F)#
Xtest <- X[idx,]
X <- X[-idx,]
Ytest <- Y[idx,]
Y <- Y[-idx,]
Ztest <- Z[idx,]
Z <- Z[-idx,]
nobs <- nrow(Y)
treat <- dplyr::recode(trtsgn, `1` = "A", `2` = "B")
treatment <- data.frame(treat[-idx])
treattest <- treat[idx]
names(treatment) <- "Treatment"
wk <- c(100, 40, 40)
#wk <- c(1, 1, 1)

#KK <- 1#numero di repliche
#res <- matrix(0, KK, 7)
#outtab <- array(0, dim = c(13, 7, KK))
#par(mfrow=c(1,1))
#for(kk in 1:KK){
modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y))
modelpriors$hP0_L0 <- diag(1, ncol(Y))
modelpriors$hP0_nu0 <- ncol(Y) + 2
modelpriors$hP0_V0 <- diag(1, ncol(Y))

alpha_DP <- 2
n_aux <- 10
vec_par <- c(0.0, 10.0, .5, 1.0, 2.0, 2.0, 0.1)
#double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
iterations <- 100
burnin <- 20
thinning <- 1

uh = F
nout <- (iterations-burnin)/thinning

time_ppm <- system.time(
  out_ppmx <- my_dm_ppmx_ct(y = Y, X = X, Xpred = Xtest,
                        z = Z, zpred = Ztest, asstreat = treatment,
                        alpha=alpha_DP, CC = n_aux, reuse = 1,
                        PPMx = 1, similarity = 1, consim=1, calibration=1,
                        similparam = vec_par, modelpriors, update_hierarchy = T,
                        iter=iterations,burn=burnin,thin=thinning, hsp = F))

print(sum(apply(round(apply(out_ppmx$isypred, c(1,2), mean))==Y, 1, sum)==3)/nobs)
#sum(apply(round(apply(out_ppmx$pi_out, c(1,2), mean))==Y, 1, sum)==3)/nobs

label <- out_ppmx$label
orig_names <- names(label)
label <- label[order(names(label))]
truelabel <- out_ppmx$asstreat
num_treat <- out_ppmx$num_treat
ntreat <- 2
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
A0
Al <- list()
for(t in 1:dim(A0)[3]){
  Al[[t]] <- A0[,,t]
}
cat("og: ", orig_names, "\n")
names(Al) <- orig_names
Al <- Al[order(names(Al))]
treatpred <- c()
for(i in 1:length(idx)){
scoreA <- sum(Al$A[i,]*wk)
scoreB <- sum(Al$B[i,]*wk)
#print(scoreA);print(scoreB);
  if(scoreA > scoreB){
    treatpred[i] <- "A"
  } else {
    treatpred[i] <- "B"
  }
}

print(treatpred)
print(treattest)

print(table(treatpred, treattest))
}
