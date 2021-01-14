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
  "C++-MN" = ran_mvnorm(theta_test, Sigma_test))
mbm
ggplot2::autoplot(mbm)

mbm <- microbenchmark::microbenchmark(
  "R-IW" = rwish(1, nu0, Sigma_test), 
  "C++-IW" = ran_wish(nu0, Sigma_test))
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
  Sn<- S0 + ( t(Y)-c(theta) )%*%t( t(Y)-c(theta) ) 
  #  Sigma<-rinvwish(1,nu0+n,solve(Sn))
  Sigma<-solve( rwish(1, nu0+n, solve(Sn)) )
  ###
  YS<-rbind(YS,rmvnorm(1,theta,Sigma)) 
  ### save results 
  THETA<-rbind(THETA,theta) ; SIGMA<-rbind(SIGMA,c(Sigma))
  #cat(s,round(theta,2),round(c(Sigma),2),"\n")
}

quantile(  SIGMA[,2]/sqrt(SIGMA[,1]*SIGMA[,4]), prob=c(.025,.5,.975) )
quantile(   THETA[,2]-THETA[,1], prob=c(.025,.5,.975) )
mean( THETA[,2]-THETA[,1])
mean( THETA[,2]>THETA[,1]) 
mean(YS[,2]>YS[,1])


##my CODE
rm(list=ls())
library(Rcpp)
source("demo/fnchoff.R")
Rcpp::sourceCpp("src/utils.cpp")

load("data/reading.RData")
Y<-reading

mu0<-c(50,50)
L0<-c( c(625,312.5,312.5,625))

nu0<-4
S0<-c( c(625,312.5,312.5,625))

n<-dim(Y)[1] ; ybar<-apply(Y,2,mean)
Sigma<-c(cov(Y)) ; THETA<-SIGMA<-NULL
YS<-NULL
#set.seed(1)
for(s in 1:5000){
  
  ###update theta
  Ln<-solve( solve(matrix(L0, 2, 2)) + n*solve(matrix(Sigma, 2, 2)) )
  mun<-Ln%*%( solve(matrix(L0, 2, 2))%*%mu0 + n*solve(matrix(Sigma, 2, 2))%*%ybar )
  theta<-c(ran_mvnorm(c(mun),c(Ln), 2)  )
  ### 
  
  ###update Sigma
  Sn<- matrix(S0, 2, 2) + ( t(Y)-c(theta) )%*%t( t(Y)-c(theta) ) 
  #  Sigma<-rinvwish(1,nu0+n,solve(Sn))
  Sigma<-solve( matrix(ran_wish(nu0+n, c(solve(Sn)), 2), 2, 2))
  #Sigma<-solve( rwish(1, nu0+n, solve(Sn)) )
  #cat("1", Sigma[1,1], "\n")
  #cat("2", Sigma[2,2], "\n")
  ###
  
  ###
  YS<-rbind(YS, c(ran_mvnorm(theta,c(Sigma), 2) ))
  ###
  
  ### save results 
  THETA<-rbind(THETA,theta) ; SIGMA<-rbind(SIGMA,c(Sigma))
  ###
  #cat(s,round(theta,2),round(c(Sigma),2),"\n")
}

quantile(  SIGMA[,2]/sqrt(SIGMA[,1]*SIGMA[,4]), prob=c(.025,.5,.975) )
quantile(   THETA[,2]-THETA[,1], prob=c(.025,.5,.975) )
mean( THETA[,2]-THETA[,1])
mean( THETA[,2]>THETA[,1]) 
mean(YS[,2]>YS[,1])

#### Figure 7.2 
par(mfrow=c(1,2),mgp=c(1.75,.75,0),mar=c(3,3,1,1))

plot.hdr2d(THETA,xlab=expression(theta[1]),ylab=expression(theta[2]) )
abline(0,1)

plot.hdr2d(YS,xlab=expression(italic(y[1])),ylab=expression(italic(y[2])), 
           xlim=c(0,100),ylim=c(0,100) )
points(Y[,1],Y[,2],pch=16,cex=.7)
abline(0,1)
