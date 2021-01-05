rm(list=ls())
library(ppmSuite)
devtools::load_all()

n <- 50
X <- cbind(rnorm(n, 0, 1), as.factor(rbinom(n, 1, 0.5)))
p <- ppmSuite::rppmx(m = n, X = X, similarity = 1, simparm = 1, M = 1)

Y <- c()
for(i in 1:p$nclus){
  for(ii in 1:n){
    if(p$Si[ii]==i){
      Y[ii] <- rnorm(1, X[i, 1]*(exp(X[i, 2])), 1)#X[which(p$Si==i),1]*log(X[which(p$Si==i), 1]) + rnorm(1, 0, .5)
    }
  }
}

X <- as.data.frame(X)
colnames(X) <- c("conx", "catx")
data <- cbind(y = as.vector(Y), X)

out <- my_ppmx(y = data$y, X = X, cohesion=1, alpha = 1.0, maug = 1, reuse = 1, 
               similarity = 2, consim = 2, calibration = 2, 
               similparam = c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1),
               modelpriors = c(0, 100^2, 0.5*sd(Y), 100), mhtune = c(1, 10),
               iter = 10000, burn = 5000, thin = 10)

pairs(data, col=out$Si[50,], pch=out$Si[50,])

par(mfrow=c(1,3))
plot(density(out$nc),type='l')
plot(density(out$mu0), type='l')
plot(density(out$sig20), type='l')

par(mfrow=c(1,2))
plot(X[,1], Y, ylab="y", xlab="x1", pch=20)
points(X[,1], apply(out$fitted,2,mean), col='blue',pch="+", cex=1.5)
acf(out$nc)
coda::effectiveSize(out$nc)