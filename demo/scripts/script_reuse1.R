rm(list=ls())
if(!is.null(dev.list())) dev.off()
# Clear console
#cat("\014") 
#devtools::load_all()
library(Rcpp)
library(rrr)
library(dplyr)
library(tidyselect)
library(mcclust)
library(mcclust.ext)
library(mclust)
library(coda)

sourceCpp("src/mvn_ppmx.cpp")
source("R/mvn_ppmx.R")
source("R/rppmx.R")

set.seed(121)

###### STUDIO 1
### Scenario a 
KK <- 2#30
res_1a <- matrix(0, KK, 8)
res_2a <- matrix(0, KK, 8)
res_3a <- matrix(0, KK, 8)
res_4a <- matrix(0, KK, 8)
res_5a <- matrix(0, KK, 8)
#set.seed(121)
K=2
myppmx <- gcd(n=100, concov = 2, K, alpha = 1)
cat("nclus: ", myppmx$nclus, "\n")
Y <- myppmx$y
colors <- c("#ebb678", "#1979a9", "#e07b39", "#69bdd2", "#80391e", "#cce7e8", 
            "#1c100b", "#042f66", "#44bcd8")
colors <- colors[myppmx$label]
if(K==2){plot(Y, pch = 16, col = colors)}
if(K==3){scatterplot3d(Y, pch = 16, color = colors, grid=TRUE, box=FALSE)}

heading <- c("nclu_post", "MSE", "lplm", "ARI", "time (sec)", "ESS", "nclu_real", "nout")
righe <- c("Aux. Cal.", " ", "Aux. Coa.", " ", "DD Cal.", " ","DD Coa", " ", "PPM", " ")
X <- myppmx$X

modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y))
modelpriors$hP0_L0 <- diag(10, ncol(Y))
modelpriors$hP0_nu0 <- nrow(Y) + 2
modelpriors$hP0_V0 <- diag(10, ncol(Y))

iterations <- 100#00
burnin <- 50#00
thinning <- 10

nout <- (iterations-burnin)/thinning


##################
## TABELLA 1R DELLO SCENARIO 1 DEL PRIMO STUDIO DI SIMULAZIONE
## n = 100, P = 2, Q = 2, K = 2, CC = 1 REUSE = TRUE
##################

for(k in 1:KK){
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 1, reuse = 1, PPMx = 1, 
                       similarity = 1, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_1a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 1, reuse = 1, PPMx = 1, 
                       similarity = 1, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_2a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 1, reuse = 1, PPMx = 1, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_3a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 1, reuse = 1, PPMx = 1, 
                       similarity = 2, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_4a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 1, reuse = 1, PPMx = 0, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_5a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
}
tab1R <- rbind(apply(res_1a, 2, mean), sqrt(apply(res_1a, 2, var)), apply(res_2a, 2, mean),
      sqrt(apply(res_2a, 2, var)), apply(res_3a, 2, mean), sqrt(apply(res_3a, 2, var)), 
      apply(res_4a, 2, mean), sqrt(apply(res_4a, 2, var)), apply(res_5a, 2, mean), 
      sqrt(apply(res_5a, 2, var)))
colnames(tab1R) <- heading
rownames(tab1R) <- righe
#save(tab1R, file = "tab1R.RData")



##################
## TABELLA 1NR DELLO SCENARIO 1 DEL PRIMO STUDIO DI SIMULAZIONE
## n = 100, P = 2, Q = 2, K = 2, CC = 1 REUSE = FALSE
##################

for(k in 1:KK){
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 1, reuse = 0, PPMx = 1, 
                       similarity = 1, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_1a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 1, reuse = 0, PPMx = 1, 
                       similarity = 1, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_2a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 1, reuse = 0, PPMx = 1, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_3a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 1, reuse = 0, PPMx = 1, 
                       similarity = 2, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_4a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 1, reuse = 0, PPMx = 0, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))

  res_5a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
}
tab1NR <- rbind(apply(res_1a, 2, mean), sqrt(apply(res_1a, 2, var)), apply(res_2a, 2, mean),
               sqrt(apply(res_2a, 2, var)), apply(res_3a, 2, mean), sqrt(apply(res_3a, 2, var)), 
               apply(res_4a, 2, mean), sqrt(apply(res_4a, 2, var)), apply(res_5a, 2, mean), 
               sqrt(apply(res_5a, 2, var)))
colnames(tab1NR) <- heading
rownames(tab1NR) <- righe
#save(tab1NR, file = "tab1NR.RData")

##################
## TABELLA 5R DELLO SCENARIO 1 DEL PRIMO STUDIO DI SIMULAZIONE
## n = 100, P = 2, Q = 2, K = 2, CC = 5 REUSE = TRUE
##################

for(k in 1:KK){
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1, 
                       similarity = 1, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_1a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1, 
                       similarity = 1, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_2a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_3a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1, 
                       similarity = 2, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_4a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 0, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_5a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
}
tab5R <- rbind(apply(res_1a, 2, mean), sqrt(apply(res_1a, 2, var)), apply(res_2a, 2, mean),
               sqrt(apply(res_2a, 2, var)), apply(res_3a, 2, mean), sqrt(apply(res_3a, 2, var)), 
               apply(res_4a, 2, mean), sqrt(apply(res_4a, 2, var)), apply(res_5a, 2, mean), 
               sqrt(apply(res_5a, 2, var)))
colnames(tab5R) <- heading
rownames(tab5R) <- righe
#save(tab5R, file = "tab5R.RData")



##################
## TABELLA 5NR DELLO SCENARIO 1 DEL PRIMO STUDIO DI SIMULAZIONE
## n = 100, P = 2, Q = 2, K = 2, CC = 1 REUSE = FALSE
##################

for(k in 1:KK){
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 0, PPMx = 1, 
                       similarity = 1, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_1a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 0, PPMx = 1, 
                       similarity = 1, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_2a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 0, PPMx = 1, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_3a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 0, PPMx = 1, 
                       similarity = 2, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_4a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 0, PPMx = 0, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_5a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
}
tab5NR <- rbind(apply(res_1a, 2, mean), sqrt(apply(res_1a, 2, var)), apply(res_2a, 2, mean),
                sqrt(apply(res_2a, 2, var)), apply(res_3a, 2, mean), sqrt(apply(res_3a, 2, var)), 
                apply(res_4a, 2, mean), sqrt(apply(res_4a, 2, var)), apply(res_5a, 2, mean), 
                sqrt(apply(res_5a, 2, var)))
colnames(tab5NR) <- heading
rownames(tab5NR) <- righe
#save(tab5NR, file = "tab5NR.RData")

##################
## TABELLA 10R DELLO SCENARIO 1 DEL PRIMO STUDIO DI SIMULAZIONE
## n = 100, P = 2, Q = 2, K = 2, CC = 10 REUSE = TRUE
##################

for(k in 1:KK){
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 10, reuse = 1, PPMx = 1, 
                       similarity = 1, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_1a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 10, reuse = 1, PPMx = 1, 
                       similarity = 1, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_2a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 10, reuse = 1, PPMx = 1, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_3a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 10, reuse = 1, PPMx = 1, 
                       similarity = 2, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_4a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 10, reuse = 1, PPMx = 0, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_5a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
}
tab10R <- rbind(apply(res_1a, 2, mean), sqrt(apply(res_1a, 2, var)), apply(res_2a, 2, mean),
               sqrt(apply(res_2a, 2, var)), apply(res_3a, 2, mean), sqrt(apply(res_3a, 2, var)), 
               apply(res_4a, 2, mean), sqrt(apply(res_4a, 2, var)), apply(res_5a, 2, mean), 
               sqrt(apply(res_5a, 2, var)))
colnames(tab10R) <- heading
rownames(tab10R) <- righe
#save(tab10R, file = "tab10R.RData")

##################
## TABELLA 10NR DELLO SCENARIO 1 DEL PRIMO STUDIO DI SIMULAZIONE
## n = 100, P = 2, Q = 2, K = 2, CC = 10 REUSE = FALSE
##################

for(k in 1:KK){
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 10, reuse = 0, PPMx = 1, 
                       similarity = 1, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_1a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 10, reuse = 0, PPMx = 1, 
                       similarity = 1, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_2a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 10, reuse = 0, PPMx = 1, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_3a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 10, reuse = 0, PPMx = 1, 
                       similarity = 2, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_4a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 10, reuse = 0, PPMx = 0, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_5a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
}
tab10NR <- rbind(apply(res_1a, 2, mean), sqrt(apply(res_1a, 2, var)), apply(res_2a, 2, mean),
                sqrt(apply(res_2a, 2, var)), apply(res_3a, 2, mean), sqrt(apply(res_3a, 2, var)), 
                apply(res_4a, 2, mean), sqrt(apply(res_4a, 2, var)), apply(res_5a, 2, mean), 
                sqrt(apply(res_5a, 2, var)))
colnames(tab10NR) <- heading
rownames(tab10NR) <- righe
#save(tab10NR, file = "tab10NR.RData")

##################
## TABELLA 30R DELLO SCENARIO 1 DEL PRIMO STUDIO DI SIMULAZIONE
## n = 100, P = 2, Q = 2, K = 2, CC = 30 REUSE = TRUE
##################

for(k in 1:KK){
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 30, reuse = 1, PPMx = 1, 
                       similarity = 1, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_1a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 30, reuse = 1, PPMx = 1, 
                       similarity = 1, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_2a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 30, reuse = 1, PPMx = 1, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_3a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 30, reuse = 1, PPMx = 1, 
                       similarity = 2, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_4a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 30, reuse = 1, PPMx = 0, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_5a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
}
tab30R <- rbind(apply(res_1a, 2, mean), sqrt(apply(res_1a, 2, var)), apply(res_2a, 2, mean),
                sqrt(apply(res_2a, 2, var)), apply(res_3a, 2, mean), sqrt(apply(res_3a, 2, var)), 
                apply(res_4a, 2, mean), sqrt(apply(res_4a, 2, var)), apply(res_5a, 2, mean), 
                sqrt(apply(res_5a, 2, var)))
colnames(tab30R) <- heading
rownames(tab30R) <- righe
#save(tab30R, file = "tab30R.RData")



##################
## TABELLA 10NR DELLO SCENARIO 1 DEL PRIMO STUDIO DI SIMULAZIONE
## n = 100, P = 2, Q = 2, K = 2, CC = 30 REUSE = FALSE
##################

for(k in 1:KK){
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 30, reuse = 0, PPMx = 1, 
                       similarity = 1, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_1a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 30, reuse = 0, PPMx = 1, 
                       similarity = 1, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_2a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 30, reuse = 0, PPMx = 1, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_3a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 30, reuse = 0, PPMx = 1, 
                       similarity = 2, consim=1, calibration=2, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_4a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
  
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 30, reuse = 0, PPMx = 0, 
                       similarity = 2, consim=1, calibration=1, 
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0), 
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), 
                       modelpriors, mhtune=c(0.5, 0.5), 
                       iter=iterations,burn=burnin,thin=thinning))
  
  res_5a[k,] <- c(unlist(postquant(lab = F)), myppmx$nclus, nout)
}
tab30NR <- rbind(apply(res_1a, 2, mean), sqrt(apply(res_1a, 2, var)), apply(res_2a, 2, mean),
                 sqrt(apply(res_2a, 2, var)), apply(res_3a, 2, mean), sqrt(apply(res_3a, 2, var)), 
                 apply(res_4a, 2, mean), sqrt(apply(res_4a, 2, var)), apply(res_5a, 2, mean), 
                 sqrt(apply(res_5a, 2, var)))
colnames(tab30NR) <- heading
rownames(tab30NR) <- righe
#save(tab30NR, file = "tab30NR.RData")

