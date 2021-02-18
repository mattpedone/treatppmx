rm(list=ls())
if(!is.null(dev.list())) dev.off()
# Clear console
#cat("\014") 
#devtools::load_all()
library(Rcpp)
library(rrr)
library(dplyr)
library(tidyselect)

sourceCpp("src/mvn_ppmx.cpp")
source("R/mvn_ppmx.R")
source("R/gendata.R")
source("R/rppmx.R")
#set.seed(24)
n=100
d <- genera_dati(n=n, P=200)

X <- d$XX

X$X1 <- as.factor(X$X1)
X$X2 <- as.factor(X$X2)

myppmx <- ran_ppmx(X=X, similarity = 2, simparm = 1, alpha=1, m0=0, s20=1,v=2, k0=10, v0=1)
myppmx

par(mfrow=c(1, 2))
possmean <- matrix(0, myppmx$nclus, 2)
for(i in 1:myppmx$nclus){
  possmean[i,] <- mvtnorm::rmvnorm(1, c(0, 0), sigma = diag(50, 2), method="svd")
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

system.time(out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 10, PPMx = 1, similarity = 2, consim=1, calibration=2, 
                   similparam=c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1, 1.0), 
                   modelpriors, mhtune=c(0.5, 0.5), 
                   iter=500,burn=250,thin=1))

#sum(myppmx$label==out$label[100,])
mean(out$nclu)

par(mfrow=c(1,2))
plot(out$nclu, type="l")
coda::effectiveSize(out$nclu)
acf(out$nclu)

table(out$label[25, ])
myppmx$nj


plot(Y, ylim=c(-40, 20), xlim=c(-20, 20))
mylab <- out$label[25,]
#par(new=T)
plot(Y[which(mylab==1), ], col = "red", ylim=c(-40, 20), xlim=c(-20, 20))
par(new=T)
plot(Y[which(mylab==2), ], col = "blue", ylim=c(-40, 20), xlim=c(-20, 20))
par(new=T)
plot(Y[which(mylab==3), ], col = "green", ylim=c(-40, 20), xlim=c(-20, 20))
par(new=T)
plot(Y[which(mylab==4), ], col = "orange", ylim=c(-40, 20), xlim=c(-20, 20))
par(new=T)
plot(Y[which(mylab==5), ], col = "black", ylim=c(-40, 20), xlim=c(-20, 20))
par(new=T)
plot(Y[which(mylab==6), ], col = "pink", ylim=c(-40, 20), xlim=c(-20, 20))

#apply(out$mu, c(1, 2), mean)
out$mu[,,25]
possmean

