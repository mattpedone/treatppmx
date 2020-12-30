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

similarityfun <- c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1)
priors <- c(0, 100^2, 0.5*sd(Y), 100)
MH <- c(1, 10)

mbm <- microbenchmark::microbenchmark("o00" = my_ppmx(y=Ytrain, X=Xtrain, Xpred=Xtest, cohesion=1, 
                                                      alpha=1.0, maug = 1, reuse = 2, similarity=2, 
                                                      consim=2, calibration=2, similparam=c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1),
                                                      modelpriors = c(0, 100^2, 0.5*sd(Y), 100),
                                                      mhtune=c(1, 10),
                                                      iter=1000,burn=500,thin=1),
                                      "o10" = my_ppmx(y=Ytrain, X=Xtrain, Xpred=Xtest, cohesion=1, 
                                                      alpha=1.0, maug = 10, reuse = 2, similarity=2, 
                                                      consim=2, calibration=2, similparam=c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1),
                                                      modelpriors = c(0, 100^2, 0.5*sd(Y), 100),
                                                      mhtune=c(1, 10),
                                                      iter=1000,burn=500,thin=1), 
                                      "o11" = my_ppmx(y=Ytrain, X=Xtrain, Xpred=Xtest, cohesion=1, 
                                                      alpha=1.0, maug = 10, reuse = 1, similarity=2, 
                                                      consim=2, calibration=2, similparam=c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1),
                                                      modelpriors = c(0, 100^2, 0.5*sd(Y), 100),
                                                      mhtune=c(1, 10),
                                                      iter=1000,burn=500,thin=1), times = 10)
mbm
ggplot2::autoplot(mbm)

out00 <- my_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, cohesion = 1, alpha = 1.0, 
                 maug = 1, reuse = 2, similarity=2, consim=2, calibration=2, 
                 similparam = similarityfun, modelpriors = priors, mhtune = MH,
                 iter = 10000, burn = 5000, thin = 1)
out10 <- my_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, cohesion = 1, alpha = 1.0, 
               maug = 10, reuse = 2, similarity=2, consim=2, calibration=2, 
               similparam = similarityfun, modelpriors = priors, mhtune = MH,
               iter = 10000, burn = 5000, thin = 1)
out11 <- my_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, cohesion = 1, alpha = 1.0, 
               maug = 10, reuse = 1, similarity=2, consim=2, calibration=2, 
               similparam = similarityfun, modelpriors = priors, mhtune = MH,
               iter = 10000, burn = 5000, thin = 1)

par(mfrow=c(1, 3))
acf(out00$nc)
acf(out10$nc)
acf(out11$nc)
coda::effectiveSize(out00$nc)
coda::effectiveSize(out10$nc)
coda::effectiveSize(out11$nc)
