#rm(list=ls())

#n = 100; P = 2; Q = 2; dim =2

genera_dati <- function(n = 100, P = 2, Q = 2, dim =2){
  discretize <- function(vec){
    vec[which(vec<0)] <- 0
    vec[which(vec>0)] <- 1
    #vec <- as.factor(vec)
    return(vec)
  }
  ncov <- Q + P
  ##prima ci stanno le binarie
  XX <- matrix(0, n, (ncov))
  
  Y <- matrix(0, n, dim)
  beta1 <- matrix(1, ncov, dim)
  beta2 <- matrix(-1.3, ncov, dim)
  beta3 <- matrix(1.3, ncov, dim)
  
  pro <- c(0.2,0.5,0.3)
  clusterlabel <- c()
  for(i in 1:n){
    u <- runif(1)
    if(u<pro[1]){
      #cat(i,"ciao1","\n")
      #X e Z li devo costruire qui
      if(Q!=0){XX[i, c(1:Q)] <- discretize(mvtnorm::rmvnorm(1, rnorm(Q, -2.1, .25)))}
      if(P!=0){XX[i, c((Q+1):ncov)] <- mvtnorm::rmvnorm(1, rnorm(P, -2.1, .25))}
      #Y[i,] <- mvtnorm::rmvnorm(1, mean=XX[i,]%*%beta1, sigma = diag(.1, nrow=dim))
      #clusterlabel[i] <- 1
    }
    else{
      if(u<(pro[1]+pro[2])){
        #cat(i,"ciao2","\n")
        if(Q!=0){XX[i, c(1:Q)] <- discretize(mvtnorm::rmvnorm(1, rnorm(Q, 0.0, .25)))}
        if(P!=0){XX[i, c((Q+1):ncov)] <- mvtnorm::rmvnorm(1, rnorm(P, 0.0, .25))}
        #Y[i,] <- mvtnorm::rmvnorm(1, mean=XX[i,]%*%beta1, sigma = diag(.1, nrow=dim))
        #clusterlabel[i] <- 2
      }
      else{
        #cat(i,"ciao3","\n")
        if(Q!=0){XX[i, c(1:Q)] <- discretize(mvtnorm::rmvnorm(1, rnorm(Q, 2.1, .25)))}
        if(P!=0){XX[i, c((Q+1):ncov)] <- mvtnorm::rmvnorm(1, rnorm(P, 2.1, .25))}
        #Y[i,] <- mvtnorm::rmvnorm(1, mean=XX[i,]%*%beta1, sigma = diag(.1, nrow=dim))
        #clusterlabel[i] <- 3
      }
    }
  }
  #XX[, seq(1:Q)] <- factor(XX[, seq(1:Q)])
  #XX <- data.frame(XX)
  #X <- XX[, -seq(1:Q)]
  return(list(XX=data.frame(XX)))
}
#
#d <- genera_dati()
#Y <- d$Y
##vedi densità bivariate
#par(mfrow=c(2, 2))
#plot(Y, xlim=c(-6, 6), ylim=c(-6, 6))
#plot(Y[which(d$cl==1),], col = "blue", xlim=c(-6, 6), ylim=c(-6, 6))
#plot(Y[which(d$cl==2),], col = "blue", xlim=c(-6, 6), ylim=c(-6, 6))
#plot(Y[which(d$cl==3),], col = "blue", xlim=c(-6, 6), ylim=c(-6, 6))

