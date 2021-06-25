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
source(file = "R/rppmx.R")
source(file = "R/dm_ppmx_ct.R")
#load("~/Dropbox/PHD/treatppmx/data/SimOutSceSimple.RData")
load("~/Dropbox/PHD/treatppmx/data/SimuOutsce2.rda")
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
k=sample(1:100, 1)
print(k)
X <- data.frame(t(mydata))[, -c(11:92)]#data.frame(mydata)#
Z <- data.frame(cbind(myx2, myx3))#data.frame(orgx)#
#Z <- apply(Z, 2, scale)
Y <- mytot[,,k]
idx <- sort(sample(1:nrow(Y), 10, replace = F))#c(23:32)#c(1:10)#
trtsgn[idx]
wk <- c(0, 40, 100)
df <- data.frame(myprob[[2]]%*%(wk)-myprob[[1]]%*%(wk))
colnames(df) <- c("Utility")
df <- cbind(Iteration = as.numeric(row.names(df)), df)
ggplot2::ggplot(df, aes(x = Iteration, y = Utility)) + geom_line() + theme_classic()

#trtsgn <- (as.numeric(myprob[[2]]%*%(wk*100)-myprob[[1]]%*%(wk*100)>0)+1)
#plot(myprob[[2]]%*%(wk*100), type = "l")
#lines(myprob[[1]]%*%(wk*100), col = "red")
#trtsgn <- (trtsgn==1)+1
treattest <- trtsgn[idx]#treat[idx]
trt <- trtsgn[-idx]
Xtest <- data.frame(X[idx,])
X <- data.frame(X[-idx,])
Ytest <- data.matrix(Y[idx,])#, nrow(Xtest), 3)
Y <- data.matrix(Y[-idx,])#, nrow(X), 3)
Ztest <- data.frame(Z[idx,])
Z <- data.frame(Z[-idx,])
nobs <- nrow(Y)

modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y))
modelpriors$hP0_L0 <- diag(10, ncol(Y))
modelpriors$hP0_nu0 <- ncol(Y) + 2
modelpriors$hP0_V0 <- diag(10, ncol(Y))

alpha_DP <- 5
n_aux <- 5
vec_par <- c(0.0, 10.0, .5, 1.0, 2.0, 2.0, 0.1)
#double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
iterations <- 1000#0
burnin <- 500
thinning <- 1#0

nout <- (iterations-burnin)/thinning
time_ppm <- system.time(
  out_ppmx <- my_dm_ppmx_ct(y = Y, X = X, Xpred = Xtest,
                        z = Z, zpred = Ztest, asstreat = trt, #treatment,
                        alpha = alpha_DP, CC = n_aux, reuse = 1,
                        PPMx = 1, similarity = 1, consim = 1, calibration = 2,
                        similparam = vec_par, modelpriors, update_hierarchy = T,
                        iter = iterations, burn = burnin, thin = thinning, hsp = T))

#traceplot for the number of clusters
df <- data.frame(t(out_ppmx$nclu))
colnames(df) <- c("t1", "t2")
df <- cbind(Index = as.numeric(row.names(df)), df)
df <- reshape2::melt(df, id.vars="Index")
#str(df)
ggplot2::ggplot(df, aes(x = Index, y = value, col = variable)) + geom_line() + theme_classic()

#a posteriori mean of prognostic covariates and some traceplots
apply(out_ppmx$beta, c(1, 2), mean)

df <- data.frame(out_ppmx$beta[1,1,])
colnames(df) <- c("b11")
df <- cbind(Iteration = as.numeric(row.names(df)), df)
ggplot2::ggplot(df, aes(x = Iteration, y = b11)) + geom_line() + theme_classic()
df <- data.frame(out_ppmx$beta[2,3,])
colnames(df) <- c("b23")
df <- cbind(Iteration = as.numeric(row.names(df)), df)
ggplot2::ggplot(df, aes(x = Iteration, y = b23)) + geom_line() + theme_classic()

#plot(eta1[2,1,], type="l")
#eta1 <- apply(eta1[,c(1:5),], c(1,2), mean)#; eta1
#eta2 <- array(unlist(out_ppmx$eta[,2]), dim = c(3, 73, nout))
#plot(eta2[3,1,], type="l")
#eta2 <- apply(eta2[,c(1:3),], c(1,2), mean)
#eta1;eta2
sum(apply(round(apply(out_ppmx$isypred, c(1,2), mean))==Y, 1, sum)==3)/nobs
#sum(apply(round(apply(out_ppmx$pi_out, c(1,2), mean))==Y, 1, sum)==3)/nobs

A0 <- out_ppmx$pipred
A0 <- array(unlist(A0), dim = c(dim(A0[[1]]), length(A0)))
A0 <- apply(A0, c(1, 2, 3), mean); A0

optrt <- as.numeric(myprob[[2]][idx,]%*%wk > myprob[[1]][idx,]%*%wk)+1
predtrt <- as.numeric(A0[,,2]%*%wk > A0[,,1]%*%wk)+1
print(optrt); print(predtrt)

#ut_sum<-sum(abs(scoreB - scoreA))
#ut_diff<- abs(as.numeric(scoreB - scoreA))

#}
