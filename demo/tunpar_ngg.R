#rm(list=ls())
#load("data/SimuOutsce2.rda")

#library(treatppmx)
library(parallel)
library(doParallel)
#source("R/countUT.R");

devtools::load_all()

K <- 10 #repliche
npat <- 152
ncolout <- (dim(mytot)[2]*2)+2+max(trtsgn)
predAPT_all<-array(0,dim=c(npat, ncolout+3,K))

wk <- c(0, 40, 100)

for(k in 1:K){
  cat("k: ", k, "\n")
  cor_all <- parallel::detectCores()-1#cores to be allocated
  registerDoParallel(cores = cor_all)

  X <- data.frame(t(mydata))[, -c(11:92)]#data.frame(mydata)#
  Z <- data.frame(cbind(myx2, myx3))
  Y <- mytot[,,k]

  modelpriors <- list()
  modelpriors$hP0_m0 <- rep(0, ncol(Y)); modelpriors$hP0_L0 <- diag(10, ncol(Y))
  modelpriors$hP0_nu0 <- ncol(Y) + 2; modelpriors$hP0_V0 <- diag(10, ncol(Y))

  n_aux <- 5 # auxiliary variable for Neal's Algorithm 8
  vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
  #double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
  iterations <- 100000;
  burnin <- 50000;
  thinning <- 10

  nout <- (iterations-burnin)/thinning
  predAPT <- c()

  myres <- foreach(sub = 1:npat, .combine = rbind) %dopar%
    {
  out_ppmx <- tryCatch(expr = my_dm_ppmx_ct(y = data.matrix(Y[-sub,]), X = data.frame(X[-sub,]), Xpred = data.frame(X[sub,]),
                            Z = data.frame(Z[-sub,]), Zpred = data.frame(Z[sub,]), asstreat = trtsgn[-sub],
                            alpha = 1, sigma = .1, CC = n_aux,
                            cohesion = 2, similarity = 2, consim = 2,
                            alphagow = 2, calibration = 2, coardegree = 2,
                            similparam = vec_par, modelpriors = modelpriors, iter = iterations,
                            burn = burnin, thin = thinning), error = function(e){FALSE})
  resvec <- c(c(apply(out_ppmx$ypred, c(1,2,3), mean)), out_ppmx$WAIC, out_ppmx$lpml, apply(out_ppmx$nclu, 1, mean))
  ifelse(is.logical(out_ppmx), return(rep(0, ncolout)), return(resvec))
    }

  ##treatment prediction with utility function ----
  A1 <- myres[,1:3]%*%wk; A2 <- myres[,4:6]%*%wk
  predAPT_all[,1,k]<-A1
  predAPT_all[,2,k]<-A2
  myt <- as.numeric(A1<A2) +1
  predAPT_all[,3,k]<-myt
  predAPT_all[,4:(ncolout+3),k]<-myres
}

my.pick <- 1:K
mywk1<-myprob[[1]]%*%wk;
mywk2<-myprob[[2]]%*%wk;
optrt<-as.numeric( mywk2> mywk1)+1;
ut.sum<-sum(abs(mywk2-mywk1));ut.diff<- abs(as.numeric(mywk2- mywk1));

#MOT
PPMXCT<-  apply(abs((predAPT_all[,3, 1:K]-optrt)),2,sum)
MOT <- c(round(mean(PPMXCT)), round(sd(PPMXCT), 1))

#MTUg
PPMXpp<-  -(2*apply(abs((predAPT_all[,3, 1:K]-optrt))*ut.diff,2,sum)-ut.sum);
MTUg <- c(round(mean(PPMXpp/ut.sum), 4), round(sd(PPMXpp/ut.sum), 4))

#NPC
PPMXCUT<-as.vector(countUT(predAPT_all));
NPC <- c(round(mean(PPMXCUT), 4), round(sd(PPMXCUT), 4))

#results
resPPMX <- rbind(MOT, MTUg,NPC)
colnames(resPPMX) <- c("mean", "sd")
resPPMX
#save(resPPMX,file="resPPMX.rda")

