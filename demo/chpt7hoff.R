#### Reading comprehension
rm(list=ls())
library(Rcpp)
library(microbenchmark)
source("demo/fnchoff.R")
Rcpp::sourceCpp("src/utils.cpp")

##COMPARISON
p <- 10
nu0 <- p + 2
theta_test <- rep(0, p)
Sigma_test <- diag(1, p, p)
mbm <- microbenchmark::microbenchmark(
  "R-MN" = rmvnorm(1, theta_test, Sigma_test), 
  "C++-MN" = ran_mvnorm(theta_test, c(Sigma_test), p))
mbm
ggplot2::autoplot(mbm)

mbm <- microbenchmark::microbenchmark(
  "R-IW" = rwish(1, nu0, Sigma_test), 
  "C++-IW" = ran_wish(nu0, c(Sigma_test), p))
mbm
ggplot2::autoplot(mbm)

##setup
#### Reading comprehension
rm(list=ls())
library(Rcpp)
source("demo/fnchoff.R")
Rcpp::sourceCpp("src/utils.cpp")

load("data/reading.RData")
Y<-reading

mu0<-c(50,50)
L0<-matrix( c(625,312.5,312.5,625),nrow=2,ncol=2)

nu0<-4
S0<-matrix( c(625,312.5,312.5,625),nrow=2,ncol=2)

n<-dim(Y)[1] ; ybar<-apply(Y,2,mean)
Sigma<-cov(Y) ; THETA<-SIGMA<-NULL
YS<-NULL
##HOFF CODE
for(s in 1:5000){
  ###update theta
  Ln<-solve( solve(L0) + n*solve(Sigma) )
  mun<-Ln%*%( solve(L0)%*%mu0 + n*solve(Sigma)%*%ybar )
  theta<-rmvnorm(1,mun,Ln)  
  ###update Sigma
  Sn<- solve(S0 + ( t(Y)-c(theta) )%*%t( t(Y)-c(theta) ) )
  Sn
  Sigma<-solve( rwish(1, nu0+n, Sn) )
  ###
  YS<-rbind(YS,rmvnorm(1,theta,Sigma)) 
  ### save results 
  THETA<-rbind(THETA,theta) ; SIGMA<-rbind(SIGMA,c(Sigma))
  }

quantile(  SIGMA[,2]/sqrt(SIGMA[,1]*SIGMA[,4]), prob=c(.025,.5,.975) )
quantile(   THETA[,2]-THETA[,1], prob=c(.025,.5,.975) )
mean( THETA[,2]-THETA[,1]); mean( THETA[,2]>THETA[,1]); mean(YS[,2]>YS[,1])


##my CODE2
rm(list=ls())
library(Rcpp)
source("demo/fnchoff.R")
Rcpp::sourceCpp("src/utils.cpp")

load("data/reading.RData")
Y<-reading

mu0<-c(50,50)
L0 <- c(matrix(c(625,312.5,312.5,625), 2, 2))
nu0<-4
S0<- c(matrix(c(625,312.5,312.5,625), 2, 2))

n<-dim(Y)[1] ; ybar<-apply(Y,2,mean)
Sigma<-c(cov(Y)) ; THETA<-SIGMA<-NULL
YS<-NULL
for(s in 1:5000){
  ###update theta
  Ln<-solve(solve(matrix(L0, 2, 2)) + n*solve(matrix(Sigma, 2, 2)))
  mun<-Ln%*%( solve(matrix(L0, 2, 2))%*%mu0 + n*solve(matrix(Sigma, 2, 2))%*%ybar )
  theta<-c(ran_mvnorm(c(mun),c(Ln), 2)  )
  ###update Sigma
  Sn<- solve(matrix(S0, 2, 2) + ( t(Y)-c(theta) )%*%t( t(Y)-c(theta)))
  #Sn
  #Sigma<-solve(matrix(ran_wish(nu0+n, c(Sn), 2, 0), 2, 2))
  Sigma<-matrix(ran_iwish(nu0+n, c(Sn), 2), 2, 2)
  ###
  YS<-rbind(YS, c(ran_mvnorm(c(theta),c(Sigma), 2) ))
  ### save results 
  THETA<-rbind(THETA,theta) ; SIGMA<-rbind(SIGMA,c(Sigma))
  }


quantile(  SIGMA[,2]/sqrt(SIGMA[,1]*SIGMA[,4]), prob=c(.025,.5,.975) )
quantile(   THETA[,2]-THETA[,1], prob=c(.025,.5,.975) )
mean( THETA[,2]-THETA[,1])
mean( THETA[,2]>THETA[,1]) 
mean(YS[,2]>YS[,1])



y <- ran_mvnorm(mu0, L0, 2)
ldmvnorm(y, mu0, matrix(L0, 2, 2))
mvtnorm::dmvnorm(c(y), mu0, matrix(L0, 2, 2), log = TRUE)
ld0 = determinant(matrix(L0, 2, 2), T)$modulus[1]
logdet(c(matrix(L0, 2, 2)), 2)
dmvnorm(y, mu0, c(matrix(L0, 2, 2)), 2, ld0, 1)

y <- ran_wish(nu0, Sig = matrix(S0, 2, 2), 2)
CholWishart::dInvWishart(matrix(y, 2, 2), nu0, matrix(S0, 2, 2), T)
LaplacesDemon::dinvwishart(matrix(y, 2, 2), nu0, matrix(S0, 2, 2), log=TRUE)
dinvwish(c(matrix(S0, 2, 2)%*%solve(matrix(y, 2, 2))), 2, 
         exp(determinant(matrix(y, 2, 2))$modulus[1]), 
         exp(determinant(matrix(S0, 2, 2))$modulus[1]), 
         nu0, 1)

#### Figure 7.2 
par(mfrow=c(1,2),mgp=c(1.75,.75,0),mar=c(3,3,1,1))

plot.hdr2d(THETA,xlab=expression(theta[1]),ylab=expression(theta[2]) )
abline(0,1)

plot.hdr2d(YS,xlab=expression(italic(y[1])),ylab=expression(italic(y[2])), 
           xlim=c(0,100),ylim=c(0,100) )
points(Y[,1],Y[,2],pch=16,cex=.7)
abline(0,1)

################################################################################
################# some checks on functions
################################################################################

rm(list=ls())
library(Rcpp)
library(geometry) 
Rcpp::sourceCpp("src/utils.cpp")

vec1 <- c(rnorm(3))
vec2 <- c(rnorm(6))
dot(vec1, vec2[c(1, 3, 5)], d = TRUE)
inner_product(vec1, 1, vec2, 2, 3)

vec1 <- matrix(rnorm(3))
norm(vec1, "2")

squared_norm(vec1, 1, 3, 0)
sqrt(squared_norm(vec1, 1, 3, 1))

rm(list=ls())
library(Rcpp)
Rcpp::sourceCpp("src/utils.cpp")

mat <- matrix(c(18, 22,  54,  42,
                22, 70,  86,  62,
                54, 86, 174, 134,
                42, 62, 134, 106), 4, 4)
#mat[lower.tri(mat)] = t(mat)[lower.tri(mat)]

chol(mat)
t(wa(c(mat), 4))
matrix(cholesky(c(mat), 4), 4, 4)

par(mfrow=c(2, 2))

mbm <- microbenchmark::microbenchmark(
  "R" = chol(mat), "work-around" = wa(c(mat), 4), "my-f" = cholesky(c(mat), 4))
mbm
ggplot2::autoplot(mbm)

determinant(mat, T)
wa_det(c(mat), 4)
logdet(c(mat), 4)

mbm <- microbenchmark::microbenchmark(
  "R" = determinant(mat, T), "work-around" = wa_det(c(mat), 4), "my-f" = logdet(c(mat), 4), times = 1000)
mbm
ggplot2::autoplot(mbm)


rm(list=ls())
library(Rcpp)
Rcpp::sourceCpp("src/utils.cpp")

mu0<-c(50,50)
L0 <- c(solve(matrix(c(625,312.5,312.5,625), 2, 2)))
nu0<-4
S0<- c(solve(matrix(c(625,312.5,312.5,625), 2, 2)))

mbm <- microbenchmark::microbenchmark(
  "prev" = ran_mvnorm_old(mu0, L0, 2), "new" = ran_mvnorm(mu0, L0, 2), times = 1000)
mbm
ggplot2::autoplot(mbm)

mbm <- microbenchmark::microbenchmark(
  "prev" = ran_iwish_old(nu0, Sig = matrix(S0, 2, 2), 2), "new" = ran_iwish(nu0, Sig = matrix(S0, 2, 2), 2), times = 1000)
mbm
ggplot2::autoplot(mbm)


rm(list=ls())
library(Rcpp)
Rcpp::sourceCpp("src/utils.cpp")

mu0<-c(50,50)
L0 <- c(solve(matrix(c(625,312.5,312.5,625), 2, 2)))

dmvnorm()