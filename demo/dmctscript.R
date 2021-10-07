  #rm(list=ls())
  #
  #library(Rcpp)
  #library(RcppArmadillo)
  #library(mvtnorm)
  #library(mcclust)
  #library(mclust)
  #library(coda)
  #library(mcclust.ext)
  ##require(ggplot2)
  ##library(reshape2)
  #
  #Rcpp::sourceCpp(file = "src/utils_ct.cpp")
  #Rcpp::sourceCpp(file = "src/dm_ppmx_ct.cpp")
  #source(file = "R/rppmx.R")
  #source(file = "R/dm_ppmx_ct.R")
  #
  #load("data/SimuOutsce2.rda")
  #data from Ma et al. Biom. J. (2017)
  #mydata contiene le covariate continue
  #myoutot contiene la risposta in ciascuna delle 100 repliche
  #versione "sintetica" di mytot
  #myprob contiene le prob di ciascun livello della risposta nei due trattamenti
  #per ciascun soggetto
  #newx sono due var normali di media 0 e var 1 e 3 (non usate??)
  #orgx contiene le due prognostiche non trasformate
  #myx2 myx3 sono le due prognostiche trasformate
  #trtsgn sono i trattamenti assegnati casualmente (trial clinico)

  devtools::load_all()

  k=sample(1:100, 1)
  print(k)
  X <- data.frame(t(mydata))[, -c(4:92)]#data.frame(mydata)#
  Z <- data.frame(cbind(myx2, myx3))#data.frame(orgx)#
  #Z <- apply(Z, 2, scale)
  Y <- mytot[,,k]
  idx <- sort(sample(1:nrow(Y), 12, replace = F))#c(35:46)#c(1:10)##
  #trtsgn[idx]
  wk <- c(0, 40, 100)
  df <- data.frame(myprob[[2]]%*%(wk)-myprob[[1]]%*%(wk))
  colnames(df) <- c("Utility")
  df <- cbind(Iteration = as.numeric(row.names(df)), df)
  ggplot2::ggplot(df, aes(x = Iteration, y = Utility)) + geom_line() + theme_classic()

  treattest <- trtsgn[idx]#treat[idx]
  trt <- trtsgn[-idx]
  Xtest <- data.frame(X[idx,])
  X <- data.frame(X[-idx,])
  Ytest <- data.matrix(Y[idx,])#, nrow(Xtest), 3)
  Y <- data.matrix(Y[-idx,])#, nrow(X), 3)
  Ztest <- data.frame(Z[idx,])
  Z <- data.frame(Z[-idx,])
  nobs <- nrow(Y)

  optrt <- as.numeric(myprob[[2]][idx,]%*%wk > myprob[[1]][idx,]%*%wk)+1; optrt

  modelpriors <- list()
  modelpriors$hP0_m0 <- rep(0, ncol(Y))
  modelpriors$hP0_L0 <- diag(10, ncol(Y))
  modelpriors$hP0_nu0 <- ncol(Y) + 2
  modelpriors$hP0_V0 <- diag(10, ncol(Y))

  alpha_DP <- 10
  n_aux <- 5
  vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
  #double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
  iterations <- 10#0000
  burnin <- 5#0000
  thinning <- 1#0

  nout <- (iterations-burnin)/thinning
  time_ppmx <- system.time(
    out_ppmx <- my_dm_ppmx_ct(y = Y, X = X, Xpred = Xtest,
                          Z = Z, Zpred = Ztest, asstreat = trt,
                          alpha = alpha_DP, sigma = .2, CC = n_aux,
                          cohesion = 2, similarity = 2, consim = 2,
                          alphagow = 2, calibration = 2, coardegree = 2,
                          similparam = vec_par, modelpriors = modelpriors, iter = iterations,
                          burn = burnin, thin = thinning))
  time_ppmx/60

  # Posterior clustering ----
  num_treat <- table(trt)
  cls <- t(as.matrix(out_ppmx$label[[1]]))[,c(1:num_treat[1])]
  psm <- comp.psm(cls)
  mc_b <- minbinder.ext(psm); max(mc_b$cl)
  mc_vi <- minVI(psm); max(mc_vi$cl)
  reord <- c()
  for(i in 1:max(mc_vi$cl)){
    reord <- c(reord, which(mc_vi$cl == i))
  }

  cls2 <- t(as.matrix(out_ppmx$label[[2]]))[,c(1:num_treat[2])]
  psm2 <- comp.psm(cls2)
  mc_b2 <- minbinder.ext(psm2); max(mc_b2$cl)
  mc_vi2 <- minVI(psm2); max(mc_vi2$cl)
  reord2 <- c()
  for(i in 1:max(mc_vi2$cl)){
    reord2 <- c(reord2, which(mc_vi2$cl == i))
  }

  # Similarity matrix ----

  melted_psm <- melt(psm)

  ggplot(data = melted_psm, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile()

  melted_psm2 <- melt(psm2)

  ggplot(data = melted_psm2, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile()

  # Co-occurence plot ----
  data <- t(out_ppmx$label[[1]])
  data <- data[,reord]
  coincidences<-sapply(1:ncol(data), function(i){ colSums(data[,i]==data) })
  ggplot(melt(coincidences), aes(Var1,Var2, fill=value)) + geom_raster() +
    scale_fill_continuous(type = "viridis")

  data <- t(out_ppmx$label[[2]])
  data <- data[,reord2]
  coincidences<-sapply(1:ncol(data), function(i){ colSums(data[,i]==data) })
  ggplot(melt(coincidences), aes(Var1,Var2, fill=value)) + geom_raster() +
    scale_fill_continuous(type = "viridis")

  # Traceplot for the number of clusters ----
  df <- data.frame(t(out_ppmx$nclu))
  colnames(df) <- c("t1", "t2")
  df <- cbind(Index = as.numeric(row.names(df)), df)
  df <- reshape2::melt(df, id.vars="Index")
  ggplot2::ggplot(df, aes(x = Index, y = value, col = variable)) + geom_line() + theme_classic()

  # A posteriori mean of prognostic covariates and some traceplots ----
  #apply(out_ppmx$beta, c(1, 2), mean)
  df <- data.frame(out_ppmx$beta[1,1,])
  colnames(df) <- c("b11")
  df <- cbind(Iteration = as.numeric(row.names(df)), df)
  ggplot2::ggplot(df, aes(x = Iteration, y = b11)) + geom_line() + theme_classic()
  df <- data.frame(out_ppmx$beta[2,3,])
  colnames(df) <- c("b23")
  df <- cbind(Iteration = as.numeric(row.names(df)), df)
  ggplot2::ggplot(df, aes(x = Iteration, y = b23)) + geom_line() + theme_classic()

  # In sample prediction (goodness-of-fit) ----
  # overall
  sum(apply(round(apply(out_ppmx$isypred, c(1,2), mean))==Y, 1, sum)==3)/nobs
  # by treatment
  sum(apply(round(apply(out_ppmx$isypred[which(trt == 1),,], c(1,2), mean))==Y[which(trt == 1),], 1, sum)==3)/sum((trt == 1))
  sum(apply(round(apply(out_ppmx$isypred[which(trt == 2),,], c(1,2), mean))==Y[which(trt == 2),], 1, sum)==3)/sum((trt == 2))

  #posterior predictive probabilities ----
  A0 <- apply(out_ppmx$ypred, c(1,2,3), mean);#A0

  #treatmente prediction with utility function ----
  optrt <- as.numeric(myprob[[2]][idx,]%*%wk > myprob[[1]][idx,]%*%wk)+1
  predtrt <- as.numeric(A0[,,2]%*%wk > A0[,,1]%*%wk)+1
  print(optrt); print(predtrt)

  sum(optrt==predtrt)/length(predtrt)
  #calcola indici paper
  #fai script x confronto

