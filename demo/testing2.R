rm(list=ls())
library(ppmSuite)
devtools::load_all()

n <- 50
X <- cbind(rnorm(n, 3, 10), as.factor(rbinom(n, 1, 0.75)))
p <- ppmSuite::rppmx(m = n, X = X, similarity = 1, simparm = 1, M = 1)
p$nclus
mv <- seq(-100, 100, length.out = p$nclus)
mv
Y <- c()
for(i in 1:p$nclus){
  for(ii in 1:n){
    if(p$Si[ii]==i){
      Y[ii] <- rnorm(1, X[i, 1]*(X[i,2]), .01)#X[which(p$Si==i),1]*log(X[which(p$Si==i), 1]) + rnorm(1, 0, .5)
    }
  }
}

#source()
#table <- read.table(file ="data/dtasim.txt", header = T)
#data <- table[sample(c(1:1000), 100),]
#Y <- data[,1]
#X <- data[,2:4]

#X[,2] <- as.factor(X[,2])
#data <- cbind(Y, X[,2], X[,1])
#pairs(data)
#set.seed(1)
train <- sample(1:length(Y), 40, replace=FALSE)
Ytrain <- as.vector(Y[train])
Ytest <- as.vector(Y[-train])
Xtrain <- data.frame(X[train,,drop=FALSE])
Xtest <- data.frame(X[-train,,drop=FALSE])


out <- my_ppmx(y=Ytrain, X=Xtrain, Xpred=Xtest, cohesion=1, 
               alpha=1.0, maug = 1, reuse = 2, similarity=1, 
               consim=1, calibration=2, similparam=c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1),
               modelpriors = c(0, 100^2, 0.5*sd(Y), 100),
               mhtune=c(1, 10),
               iter=100,burn=50,thin=1)


#plot(density(out$mu[,1:5000]),type='l')
#plot(density(out$sig2[,1:10]),type='l')
plot(density(out$nc),type='l')
plot(density(out$mu0), type='l')
plot(density(out$sig20), type='l')

plot(Xtrain[,1], Ytrain, ylab="y", xlab="x1", pch=20)
points(Xtrain[,1], apply(out$fitted,2,mean), col='blue',pch="+", cex=1.5)
plot(Xtest[,1], Ytest, ylab="y", xlab="x1",pch=20)
points(Xtest[,1], apply(out$ppred,2,mean), col='blue',pch="+", cex=1.5)

