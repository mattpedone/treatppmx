devtools::load_all()
library(ppmSuite) 
data(bear)
# plot length, sex, and weight of bears
ck <- c(4,3,2)
pairs(bear[,ck])
# response is length
Y <- bear$weight
# Continuous Covariate is chest
# Categorical covariate is sex
X <- bear[,c("length", "sex")]
X$sex <- as.factor(X$sex)

# Randomly partition data into 44 training and 10 testing
set.seed(1)
trainObs <- sample(1:length(Y),44, replace=FALSE)
Ytrain <- Y[trainObs]
Ytest <- Y[-trainObs]
Xtrain <- X[trainObs,,drop=FALSE]
Xtest <- X[-trainObs,,drop=FALSE]


out <- my_ppmx(y=Ytrain, X=Xtrain, Xpred=Xtest, cohesion=1, alpha=1.0, maug = 3,
        similarity=2, consim=2, calibration=2,
        similparam=c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1),
        modelpriors = c(0, 100^2, 0.5*sd(Y), 100),
        mhtune=c(1, 10),
        iter=100000,burn=50000,thin=50)
pairs(bear[trainObs,ck], col=out$Si[1000,], pch=out$Si[1000,])

## plot MCMC iterats
plot(density(out$mu[,1:10]),type='l')
plot(density(out$sig2[,1:10]),type='l')
plot(density(out$nc),type='l')
plot(density(out$mu0), type='l')
plot(density(out$sig20), type='l')

plot(Xtrain[,1], Ytrain, ylab="weight", xlab="length", pch=20)
points(Xtrain[,1], apply(out$fitted,2,mean), col='blue',pch="+", cex=1.5)
#points(Xtrain[,1], apply(out2$fitted,2,mean), col='red',pch=2, cex=1)
#legend(x="topleft",legend=c("Observed","PPMx","PPM"), col=c("black","red","blue"),pch=c(20,3,2))
plot(Xtest[,1], Ytest, ylab="weight", xlab="length",pch=20)
points(Xtest[,1], apply(out$ppred,2,mean), col='blue',pch="+", cex=1.5)
#points(Xtest[,1], apply(out2$ppred,2,mean), col='red',pch=2, cex=1)
#legend(x="topleft",legend=c("Observed","PPMx","PPM"), col=c("black","blue","red"),pch=c(20,3,2))
#par(oldpar)
#
#
