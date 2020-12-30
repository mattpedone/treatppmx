rm(list=ls())
devtools::load_all()

#costruisci dati seguendo modello generativo

set.seed(1)
train <- sample(1:length(Y), 75, replace=FALSE)
Ytrain <- Y[train]
Ytest <- Y[-train]
Xtrain <- X[train,,drop=FALSE]
Xtest <- X[-train,,drop=FALSE]

out <- my_ppmx(y=Ytrain, X=Xtrain, Xpred=Xtest, cohesion=1, 
               alpha=1.0, maug = 1, reuse = 2, similarity=2, 
               consim=1, calibration=2, similparam=c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1),
               modelpriors = c(0, 100^2, 0.5*sd(Y), 100),
               mhtune=c(1, 10),
               iter=10000,burn=5000,thin=1)


plot(density(out$mu[,1:10]),type='l')
plot(density(out$sig2[,1:10]),type='l')
plot(density(out$nc),type='l')
plot(density(out$mu0), type='l')
plot(density(out$sig20), type='l')

plot(Xtrain[,1], Ytrain, ylab="y", xlab="x1", pch=20)
points(Xtrain[,1], apply(out$fitted,2,mean), col='blue',pch="+", cex=1.5)
plot(Xtest[,1], Ytest, ylab="y", xlab="x1",pch=20)
points(Xtest[,1], apply(out$ppred,2,mean), col='blue',pch="+", cex=1.5)

