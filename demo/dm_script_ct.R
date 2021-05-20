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
library(dplyr)
library(lessR)

Rcpp::sourceCpp(file = "src/utils_ct.cpp")
Rcpp::sourceCpp(file = "src/dm_ppmx_ct.cpp")
source(file = "R/dm_ppmx_ct.R")
source(file = "R/rppmx.R")
set.seed(121)

KK <- 1#numero di repliche
res <- matrix(0, KK, 7)
outtab <- array(0, dim = c(13, 7, KK))
par(mfrow=c(1,1))
#for(kk in 1:KK){
  mydata <- scenario2_ct()
  Y <- mydata$y
  X <- mydata$X
  Ytrain <- mydata$Ytrain
  Ytest <- mydata$Ytest
  Xtrain <- mydata$Xtrain
  Xtest <- mydata$Xtest
  Z <- mydata$Z
  Ztrain <- mydata$Ztrain
  Ztest <- mydata$Ztest
  asstreat <- data.frame(Treatment = mydata$treat)

  for(j in 1:mydata$nclus){
    cat(apply(matrix(Y[which(mydata$labeltrain == j),], ncol=ncol(Y)), 2, sum), "\n")
  }

  modelpriors <- list()
  modelpriors$hP0_m0 <- rep(0, ncol(Y))
  modelpriors$hP0_L0 <- diag(1, ncol(Y))
  modelpriors$hP0_nu0 <- ncol(Y) + 2
  modelpriors$hP0_V0 <- diag(1, ncol(Y))

  alpha_DP <- 2
  n_aux <- 5
  vec_par <- c(0.0, 10.0, .5, 1.0, 2.0, 2.0, 0.1)
  #double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
  iterations <- 520
  burnin <- 20
  thinning <- 1

  uh = T
  nout <- (iterations-burnin)/thinning

  # PPM
  time_ppm <- system.time(
    out_ppm <- my_dm_ppmx_ct(y = Ytrain, X = Xtrain, Xpred = Xtest,
                          z = Ztrain, zpred = Ztest, asstreat = asstreat, alpha=2,
                          CC = n_aux, reuse = 1,
                          PPMx = 0, similarity = 1, consim=1, calibration=0,
                          similparam = vec_par, modelpriors, update_hierarchy = F,
                          iter=iterations,burn=burnin,thin=thinning))
  ppm <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppm, minbinder = T,
                      data = mydata, plot = F)

  # PPMx No Calibration Auxiliary Similarity NN
  time_ppmx0_aux <- system.time(
    out_ppmx0_aux <- my_dm_ppmx_ct(y = Ytrain, X = Xtrain, Xpred = Xtest,
                                z = Ztrain, zpred = Ztest, asstreat = traintreat, alpha=alpha_DP,
                                CC = n_aux, reuse = 1,
                                PPMx = 1, similarity = 1, consim=1, calibration=0,
                                similparam = vec_par, modelpriors,  update_hierarchy = uh,
                                iter=iterations,burn=burnin,thin=thinning))
  ppmx0_aux <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx0_aux,
                            data = mydata, plot = F)

  probs_ppm <- apply(out_ppm$pipred, c(1, 2), mean)
  probs_ppmx <- apply(out_ppmx0_aux$pipred, c(1, 2), mean)

  colnames(Ytest) <- c("a_true", "b_true", "c_true", "d_true")
  colnames(probs_ppm) <- c("a_pred_ppm", "b_pred_ppm", "c_pred_ppm", "d_pred_ppm")
  colnames(probs_ppmx) <- c("a_pred_ppmx", "b_pred_ppmx", "c_pred_ppmx", "d_pred_ppmx")
  final_df <- data.frame(cbind(Ytest, probs_ppm, probs_ppmx))

  roc_res <- multi_roc(final_df, force_diag = T)
  plot_roc_df <- plot_roc_data(roc_res)

  ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
    geom_path(aes(color = Group, linetype=Method), size=0.5)

  out_ppmx0_aux$acc_rate_beta/iterations

  # PPMx No Calibration Auxiliary Similarity NNIG
  time_ppmx0_aux_IG <- system.time(
    out_ppmx0_aux_IG <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest,
                                   z = Ztrain, zpred = Ztest, alpha=alpha_DP,
                                   CC = n_aux, reuse = 1,
                                   PPMx = 1, similarity = 1, consim=2, calibration=0,
                                   similparam = vec_par, modelpriors,  update_hierarchy = uh,
                                   iter=iterations,burn=burnin,thin=thinning))
  ppmx0_aux_IG <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx0_aux_IG,
                               data = mydata, plot = F)

  # PPMx No Calibration Double Dipper similarity
  time_ppmx0_dd <- system.time(
    out_ppmx0_dd <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest,
                               z = Ztrain, zpred = Ztest, alpha=alpha_DP,
                               CC = n_aux, reuse = 1,
                               PPMx = 1, similarity = 2, consim=1, calibration=0,
                               similparam = vec_par, modelpriors, update_hierarchy = uh,
                               iter=iterations,burn=burnin,thin=thinning))
  ppmx0_dd <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx0_dd,
                           data = mydata, plot = F)

  # PPMx No Calibration Double Dipper similarity NNIG
  time_ppmx0_dd_IG <- system.time(
    out_ppmx0_dd_IG <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest,
                                  z = Ztrain, zpred = Ztest, alpha=alpha_DP,
                                  CC = n_aux, reuse = 1,
                                  PPMx = 1, similarity = 2, consim=2, calibration=0,
                                  similparam = vec_par, modelpriors, update_hierarchy = uh,
                                  iter=iterations,burn=burnin,thin=thinning))
  ppmx0_dd_IG <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx0_dd_IG,
                              data = mydata, plot = F)

  # PPMx Calibrated Auxiliary similarity
  time_ppmx1_aux <- system.time(
    out_ppmx1_aux <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest,
                                z = Ztrain, zpred = Ztest, alpha=alpha_DP,
                                CC = n_aux, reuse = 1,
                                PPMx = 1, similarity = 1, consim=1, calibration=1,
                                similparam = vec_par, modelpriors, update_hierarchy = uh,
                                iter=iterations,burn=burnin,thin=thinning))
  ppmx1_aux <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx1_aux,
                            data = mydata, plot = F)

  # PPMx Calibrated Auxiliary similarity NNIG
  time_ppmx1_aux_IG <- system.time(
    out_ppmx1_aux_IG <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest,
                                   z = Ztrain, zpred = Ztest, alpha=alpha_DP,
                                   CC = n_aux, reuse = 1,
                                   PPMx = 1, similarity = 1, consim=2, calibration=1,
                                   similparam = vec_par, modelpriors, update_hierarchy = uh,
                                   iter=iterations,burn=burnin,thin=thinning))
  ppmx1_aux_IG <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx1_aux_IG,
                               data = mydata, plot = F)

  # PPMx Calibrated Double Dipper similarity
  time_ppmx1_dd <- system.time(
    out_ppmx1_dd <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest,
                               z = Ztrain, zpred = Ztest, alpha=alpha_DP,
                               CC = n_aux, reuse = 1,
                               PPMx = 1, similarity = 2, consim=1, calibration=1,
                               similparam = vec_par, modelpriors, update_hierarchy = uh,
                               iter=iterations,burn=burnin,thin=thinning))
  ppmx1_dd <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx1_dd,
                           data = mydata, plot = F)

  # PPMx Calibrated Double Dipper similarity NNIG
  time_ppmx1_dd_IG <- system.time(
    out_ppmx1_dd_IG <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest,
                                  z = Ztrain, zpred = Ztest, alpha=alpha_DP,
                                  CC = n_aux, reuse = 1,
                                  PPMx = 1, similarity = 2, consim=2, calibration=1,
                                  similparam = vec_par, modelpriors, update_hierarchy = uh,
                                  iter=iterations,burn=burnin,thin=thinning))
  ppmx1_dd_IG <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx1_dd_IG,
                              data = mydata, plot = F)

  # PPMx Coarsened Auxiliary similarity
  time_ppmx2_aux <- system.time(
    out_ppmx2_aux <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest,
                                z = Ztrain, zpred = Ztest, alpha=alpha_DP,
                                CC = n_aux, reuse = 1,
                                PPMx = 1, similarity = 1, consim=1, calibration=2,
                                similparam = vec_par, modelpriors, update_hierarchy = uh,
                                iter=iterations,burn=burnin,thin=thinning))
  ppmx2_aux <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx2_aux,
                            data = mydata, plot = F)

  # PPMx Coarsened Auxiliary similarity NNIG
  time_ppmx2_aux_IG <- system.time(
    out_ppmx2_aux_IG <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest,
                                   z = Ztrain, zpred = Ztest, alpha=alpha_DP,
                                   CC = n_aux, reuse = 1,
                                   PPMx = 1, similarity = 1, consim=2, calibration=2,
                                   similparam = vec_par, modelpriors, update_hierarchy = uh,
                                   iter=iterations,burn=burnin,thin=thinning))
  ppmx2_aux_IG <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx2_aux_IG,
                               data = mydata, plot = F)

  # PPMx Coarsened Double Dipper similarity
  time_ppmx2_dd <- system.time(
    out_ppmx2_dd <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest,
                               z = Ztrain, zpred = Ztest, alpha=alpha_DP,
                               CC = n_aux, reuse = 1,
                               PPMx = 1, similarity = 2, consim=1, calibration=2,
                               similparam = vec_par, modelpriors, update_hierarchy = uh,
                               iter=iterations,burn=burnin,thin=thinning))
  ppmx2_dd <- postquant_dm(y = Ytrain, yp = Ytest, output = out_ppmx2_dd,
                           data = mydata, plot = F)

  # PPMx Coarsened Double Dipper similarity NNIG
  time_ppmx2_dd_IG <- system.time(
    out_ppmx2_dd_IG <- my_dm_ppmx(y = Ytrain, X = Xtrain, Xpred = Xtest,
                                  z = Ztrain, zpred = Ztest, alpha=alpha_DP,
                                  CC = n_aux, reuse = 1,
                                  PPMx = 1, similarity = 2, consim=2, calibration=2,
                                  similparam = vec_par, modelpriors, update_hierarchy = uh,
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
  kk=1
  outtab[,,kk] <- matrix(data.matrix(unlist(tab)), ncol=7)
#}
#sum(unlist(tab[,7]))/60

bestauc <- which(tab[,6] == max(unlist(tab[,6])))
ppmxs <- rownames(tab)
plot_auc(out_ppm, eval(parse(text =paste0("out_", ppmxs[bestauc]))))


mean_rep <- apply(outtab, c(1, 2), mean)
colnames(mean_rep) <- colnames(tab)
rownames(mean_rep) <- rownames(tab)
mean_rep

sd_rep <- apply(outtab, c(1, 2), var)
colnames(sd_rep) <- colnames(tab)
rownames(sd_rep) <- rownames(tab)
sqrt(sd_rep)
