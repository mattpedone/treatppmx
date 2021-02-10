devtools::load_all()
library(ppmSuite) 
data(bear)

ck <- c(4,3,2)
pairs(bear[,ck])
Y <- bear$weight
X <- bear[,c("length", "sex")]
X$sex <- as.factor(X$sex)

out <- my_ppmx(y = Y, X = X, cohesion=1, 
               alpha=1.0, maug = 1, #reuse = 1, 
               similarity=2, 
               consim=2, calibration=2, similparam=c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1),
                modelpriors = c(0, 100^2, 0.5*sd(Y), 100),
                mhtune=c(1, 10),
                iter=50000,burn=25000,thin=50)

par(mfrow=c(1,2))
plot(out$nc, type="l")
coda::effectiveSize(out$nc)
acf(out$nc)

pairs(bear[,ck], col=out$Si[50,], pch=out$Si[50,])

## plot MCMC iterats
par(mfrow=c(1,3))
plot(density(out$mu[,1:dim(out$nc)[1]]),type='l')
plot(density(out$sig2[,1:dim(out$nc)[1]]),type='l')
plot(density(out$nc),type='l')
par(mfrow=c(1,3))
plot(density(out$mu0), type='l')
plot(density(out$sig20), type='l')
plot(X[,1], Y, ylab="weight", xlab="length", pch=20)
points(X[,1], apply(out$fitted,2,mean), col='blue',pch="+", cex=1.5)

rm(list=ls())
devtools::load_all()
library(rrr)
library(dplyr)
library(tidyselect)
data(tobacco)
tobacco <- as_data_frame(tobacco)
tobacco_x <- tobacco[, -c(1:3)]

tobacco_y <- tobacco[, c(1:3)]

#x <- as.matrix(tobacco_x)
#y <- as.matrix(tobacco_y)

modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, 3)
modelpriors$hP0_L0 <- diag(1, 3)
modelpriors$hP0_nu0 <- nrow(tobacco_y) + 2
modelpriors$hP0_V0 <- diag(1, 3)

mr <- my_mvn_ppmx(y = tobacco_y, X = tobacco_x, alpha=1, CC = 1, similarity = 2, consim=1, calibration=2,
                        similparam=c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1, 1.0),
                        modelpriors,
                        mhtune=c(0.5, 0.5),
                        iter=100,burn=50,thin=10)

mr$nclu
##è come se ci fosse un disallineamento tra curr_clu e nj_curr mi sembra che quando dei soggetti da un cluster vengono spostati
## il cluster iniziale non sia svuotato. quindi il conteggio dei cluster viene correttamente "scalato"
##la loro cardinalità non è aggiornata però.


