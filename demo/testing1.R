devtools::load_all()
library(ppmSuite) 
data(bear)

ck <- c(4,3,2)
pairs(bear[,ck])
Y <- bear$weight
X <- bear[,c("length", "sex")]
X$sex <- as.factor(X$sex)

similarityfun <- c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1)
priors <- c(0, 100^2, 0.5*sd(Y), 100)
MH <- c(1, 10)

mbm <- microbenchmark::microbenchmark(
  "o10" = my_ppmx(y = Y, X = X, cohesion = 1, alpha = 1.0, maug = 1, reuse = 2, 
                  similarity = 2, consim = 2, calibration = 2, 
                  similparam = similarityfun, modelpriors = priors, mhtune = MH, 
                  iter = 1000,burn = 500,thin = 1), 
  "o11" = my_ppmx(y = Y, X = X, cohesion = 1, alpha = 1.0, maug = 1, reuse = 1, 
                  similarity = 2, consim = 2, calibration = 2, 
                  similparam = similarityfun, modelpriors = priors, mhtune = MH, 
                  iter = 1000,burn = 500,thin = 1), 
  "o100" = my_ppmx(y = Y, X = X, cohesion = 1, alpha = 1.0, maug = 10, reuse = 2, 
                  similarity = 2, consim = 2, calibration = 2, 
                  similparam = similarityfun, modelpriors = priors, mhtune = MH, 
                  iter = 1000,burn = 500,thin = 1), 
  "o101" = my_ppmx(y = Y, X = X, cohesion = 1, alpha = 1.0, maug = 10, reuse = 1, 
                  similarity = 2, consim = 2, calibration = 2, 
                  similparam = similarityfun, modelpriors = priors, mhtune = MH, 
                  iter = 1000,burn = 500,thin = 1), times = 10)
mbm
ggplot2::autoplot(mbm)

out10 <- my_ppmx(y = Y, X = X, cohesion = 1, alpha = 1.0, maug = 1, reuse = 2, 
                 similarity = 2, consim = 2, calibration = 2, 
                 similparam = similarityfun, modelpriors = priors, mhtune = MH, 
                 iter = 10000,burn = 5000,thin = 1)

out11 <- my_ppmx(y = Y, X = X, cohesion = 1, alpha = 1.0, maug = 1, reuse = 1, 
                 similarity = 2, consim = 2, calibration = 2, 
                 similparam = similarityfun, modelpriors = priors, mhtune = MH, 
                 iter = 10000,burn = 5000,thin = 1)

out100 <- my_ppmx(y = Y, X = X, cohesion = 1, alpha = 1.0, maug = 10, reuse = 2, 
                 similarity = 2, consim = 2, calibration = 2, 
                 similparam = similarityfun, modelpriors = priors, mhtune = MH, 
                 iter = 10000,burn = 5000,thin = 1)

out101 <- my_ppmx(y = Y, X = X, cohesion = 1, alpha = 1.0, maug = 10, reuse = 1, 
                 similarity = 2, consim = 2, calibration = 2, 
                 similparam = similarityfun, modelpriors = priors, mhtune = MH, 
                 iter = 10000,burn = 5000,thin = 1)

par(mfrow=c(1, 2))
acf(out10$nc)
acf(out11$nc)
acf(out100$nc)
acf(out101$nc)
coda::effectiveSize(out10$nc)
coda::effectiveSize(out11$nc)
coda::effectiveSize(out100$nc)
coda::effectiveSize(out101$nc)
