rm(list=ls())

#devtools::load_all()
library(mvtnorm)
library(mcclust)
library(mclust)
library(coda)
library(mcclust.ext)
library(pROC)
library(multiROC)
require(ggplot2)

Rcpp::sourceCpp(file = "src/utils.cpp")
Rcpp::sourceCpp(file = "src/dm_ppmx.cpp")
source(file = "R/dm_ppmx.R")
source(file = "R/rppmx.R")

KK <- 1#numero di repliche
res <- matrix(0, KK, 7)
par(mfrow=c(1,1))

mydata <- scenario2()
Y <- mydata$y
X <- mydata$X
Ytrain <- mydata$Ytrain
Ytest <- mydata$Ytest
Xtrain <- mydata$Xtrain
Xtest <- mydata$Xtest

for(j in 1:mydata$nclus){
  cat(apply(matrix(Y[which(mydata$labeltrain == j),], ncol=ncol(Y)), 2, sum), "\n")
}

modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y))
modelpriors$hP0_L0 <- diag(1, ncol(Y))
modelpriors$hP0_nu0 <- ncol(Y) + 2
modelpriors$hP0_V0 <- diag(1, ncol(Y))

alpha_DP <- 1
n_aux <- 5
vec_par <- c(0.0, 10.0, .5, 1.0, 2.0, 2.0, 0.1)
#double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
iterations <- 52000
burnin <- 2000
thinning <- 10

nout <- (iterations-burnin)/thinning

# PPM
time_ppm <- system.time(
  out_ppm <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, alpha=alpha_DP,
                        CC = n_aux, reuse = 1,
                        PPMx = 0, similarity = 1, consim=1, calibration=0,
                        similparam = vec_par, modelpriors, update_hierarchy = F,
                        iter=iterations,burn=burnin,thin=thinning))
ppm <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppm,
                    data = mydata, plot = F)

# PPMx No Calibration Auxiliary Similarity NN
time_ppmx0_aux <- system.time(
  out_ppmx0_aux <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, alpha=alpha_DP,
                              CC = n_aux, reuse = 1,
                              PPMx = 1, similarity = 1, consim=1, calibration=0,
                              similparam = vec_par, modelpriors,  update_hierarchy = F,
                              iter=iterations,burn=burnin,thin=thinning))
ppmx0_aux <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx0_aux,
                          data = mydata, plot = F)

# PPMx No Calibration Auxiliary Similarity NNIG
time_ppmx0_aux_IG <- system.time(
  out_ppmx0_aux_IG <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, alpha=alpha_DP,
                                 CC = n_aux, reuse = 1,
                              PPMx = 1, similarity = 1, consim=2, calibration=0,
                              similparam = vec_par, modelpriors,  update_hierarchy = F,
                              iter=iterations,burn=burnin,thin=thinning))
ppmx0_aux_IG <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx0_aux_IG,
                             data = mydata, plot = F)

# PPMx No Calibration Double Dipper similarity
time_ppmx0_dd <- system.time(
  out_ppmx0_dd <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, alpha=alpha_DP,
                             CC = n_aux, reuse = 1,
                             PPMx = 1, similarity = 2, consim=1, calibration=0,
                             similparam = vec_par, modelpriors, update_hierarchy = F,
                             iter=iterations,burn=burnin,thin=thinning))
ppmx0_dd <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx0_dd,
                         data = mydata, plot = F)

# PPMx No Calibration Double Dipper similarity NNIG
time_ppmx0_dd_IG <- system.time(
  out_ppmx0_dd_IG <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, alpha=alpha_DP,
                                CC = n_aux, reuse = 1,
                             PPMx = 1, similarity = 2, consim=2, calibration=0,
                             similparam = vec_par, modelpriors, update_hierarchy = F,
                             iter=iterations,burn=burnin,thin=thinning))
ppmx0_dd_IG <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx0_dd_IG,
                            data = mydata, plot = F)

# PPMx Calibrated Auxiliary similarity
time_ppmx1_aux <- system.time(
  out_ppmx1_aux <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, alpha=alpha_DP,
                              CC = n_aux, reuse = 1,
                              PPMx = 1, similarity = 1, consim=1, calibration=1,
                              similparam = vec_par, modelpriors, update_hierarchy = F,
                              iter=iterations,burn=burnin,thin=thinning))
ppmx1_aux <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx1_aux,
                          data = mydata, plot = F)

# PPMx Calibrated Auxiliary similarity NNIG
time_ppmx1_aux_IG <- system.time(
  out_ppmx1_aux_IG <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, alpha=alpha_DP,
                                 CC = n_aux, reuse = 1,
                              PPMx = 1, similarity = 1, consim=2, calibration=1,
                              similparam = vec_par, modelpriors, update_hierarchy = F,
                              iter=iterations,burn=burnin,thin=thinning))
ppmx1_aux_IG <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx1_aux_IG,
                             data = mydata, plot = F)

# PPMx Calibrated Double Dipper similarity
time_ppmx1_dd <- system.time(
  out_ppmx1_dd <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, alpha=alpha_DP,
                             CC = n_aux, reuse = 1,
                             PPMx = 1, similarity = 2, consim=1, calibration=1,
                             similparam = vec_par, modelpriors, update_hierarchy = F,
                             iter=iterations,burn=burnin,thin=thinning))
ppmx1_dd <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx1_dd,
                         data = mydata, plot = F)

# PPMx Calibrated Double Dipper similarity NNIG
time_ppmx1_dd_IG <- system.time(
  out_ppmx1_dd_IG <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, alpha=alpha_DP,
                                CC = n_aux, reuse = 1,
                             PPMx = 1, similarity = 2, consim=2, calibration=1,
                             similparam = vec_par, modelpriors, update_hierarchy = F,
                             iter=iterations,burn=burnin,thin=thinning))
ppmx1_dd_IG <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx1_dd_IG,
                            data = mydata, plot = F)

# PPMx Coarsened Auxiliary similarity
time_ppmx2_aux <- system.time(
  out_ppmx2_aux <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, alpha=alpha_DP,
                              CC = n_aux, reuse = 1,
                              PPMx = 1, similarity = 1, consim=1, calibration=2,
                              similparam = vec_par, modelpriors, update_hierarchy = F,
                              iter=iterations,burn=burnin,thin=thinning))
ppmx2_aux <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx2_aux,
                          data = mydata, plot = F)

# PPMx Coarsened Auxiliary similarity NNIG
time_ppmx2_aux_IG <- system.time(
  out_ppmx2_aux_IG <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, alpha=alpha_DP,
                                 CC = n_aux, reuse = 1,
                              PPMx = 1, similarity = 1, consim=2, calibration=2,
                              similparam = vec_par, modelpriors, update_hierarchy = F,
                              iter=iterations,burn=burnin,thin=thinning))
ppmx2_aux_IG <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx2_aux_IG,
                             data = mydata, plot = F)

# PPMx Coarsened Double Dipper similarity
time_ppmx2_dd <- system.time(
  out_ppmx2_dd <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, alpha=alpha_DP,
                             CC = n_aux, reuse = 1,
                             PPMx = 1, similarity = 2, consim=1, calibration=2,
                             similparam = vec_par, modelpriors, update_hierarchy = F,
                             iter=iterations,burn=burnin,thin=thinning))
ppmx2_dd <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx2_dd,
                         data = mydata, plot = F)

# PPMx Coarsened Double Dipper similarity NNIG
time_ppmx2_dd_IG <- system.time(
  out_ppmx2_dd_IG <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest, alpha=alpha_DP,
                                CC = n_aux, reuse = 1,
                             PPMx = 1, similarity = 2, consim=2, calibration=2,
                             similparam = vec_par, modelpriors, update_hierarchy = F,
                             iter=iterations,burn=burnin,thin=thinning))
ppmx2_dd_IG <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx2_dd_IG,
                            data = mydata, plot = F)

tab <- rbind(ppm, ppmx0_aux, ppmx0_aux_IG, ppmx0_dd, ppmx0_dd_IG, ppmx1_aux,
             ppmx1_aux_IG, ppmx1_dd, ppmx1_dd_IG, ppmx2_aux, ppmx2_aux_IG,
             ppmx2_dd, ppmx2_dd_IG)

time <- c(time_ppm[3], time_ppmx0_aux[3], time_ppmx0_aux_IG[3], time_ppmx0_dd[3],
          time_ppmx0_dd_IG[3], time_ppmx1_aux[3], time_ppmx1_aux_IG[3],
          time_ppmx1_dd[3], time_ppmx1_dd_IG[3], time_ppmx2_aux[3], time_ppmx2_aux_IG[3],
          time_ppmx2_dd[3], time_ppmx2_dd_IG[3])
tab <- cbind(tab[,-4], time)
tab
sum(unlist(tab[,7]))/60

bestauc <- which(tab[,6] == max(unlist(tab[,6])))
ppmxs <- rownames(tab)
plot_auc(out_ppm, eval(parse(text =paste0("out_", ppmxs[bestauc]))))

#postquant_dm(y = Y, output = out_ppm, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx0_aux, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx0_aux_IG, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx0_dd, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx0_dd_IG, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx1_aux, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx1_aux_IG, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx1_dd, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx1_dd_IG, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx2_aux, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx2_aux_IG, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx2_dd, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx2_dd_IG, data = mydata, plot = T)
