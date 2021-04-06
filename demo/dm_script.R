rm(list=ls())

devtools::load_all()

KK <- 1#numero di repliche
res <- matrix(0, KK, 7)
par(mfrow=c(1,1))
K=4#dimensioni
mydata <- gcd_dm(n_obs = 100, concov = 2, K, similarity = 1, simparm = 1,
              alpha = 1, m0 = 0, s20 = 1, v = 2, k0 = 10, v0 = 1, plot = F)
mydata$possmean
cat("nclus: ", mydata$nclus, "\n")
#mydata <- gcd_dm_simpl(100)
Y <- mydata$y

for(j in 1:mydata$nclus){
  cat(apply(matrix(Y[which(mydata$label == j),], ncol=K), 2, sum), "\n")
}

heading <- c("nclu_post", "MSE", "lplm", "ARI", "ESS", "nclu_real", "nout")
#righe <- c("Reuse", " ", "No ", " ", "DD Cal.", " ","DD Coa", " ", "PPM", " ")
X <- mydata$X
#plot(X[,1], X[, 2])
modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y))
modelpriors$hP0_L0 <- diag(1, ncol(Y))
modelpriors$hP0_nu0 <- nrow(Y) + 20
modelpriors$hP0_V0 <- diag(100, ncol(Y))

iterations <- 100
burnin <- 0
thinning <- 1

nout <- (iterations-burnin)/thinning

mytime_ppm <- system.time(
  out_ppm <- my_dm_ppmx(y = Y, X = X, alpha=1, CC = 1, reuse = 0, PPMx = 0, similarity = 1, consim=1, calibration=0,
                     similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), modelpriors, mhtune=c(0.5, 0.5),
                     iter=iterations,burn=burnin,thin=thinning))
po_ppm <- postquant_dm(y = Y, output = out_ppm, data = mydata, plot = F)

mytime_ppmx0_aux <- system.time(
  out_ppmx0_aux <- my_dm_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1, similarity = 1, consim=1, calibration=0,
                    similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), modelpriors, mhtune=c(0.5, 0.5),
                    iter=iterations,burn=burnin,thin=thinning))
po_ppmx0_aux <- postquant_dm(y = Y, output = out_ppmx0_aux, data = mydata, plot = F)

mytime_ppmx0_dd <- system.time(
  out_ppmx0_dd <- my_dm_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1, similarity = 2, consim=1, calibration=0,
                              similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), modelpriors, mhtune=c(0.5, 0.5),
                              iter=iterations,burn=burnin,thin=thinning))
po_ppmx0_dd <- postquant_dm(y = Y, output = out_ppmx0_dd, data = mydata, plot = F)

mytime_ppmx1_aux <- system.time(
  out_ppmx1_aux <- my_dm_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1, similarity = 1, consim=1, calibration=1,
                    similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), modelpriors, mhtune=c(0.5, 0.5),
                    iter=iterations,burn=burnin,thin=thinning))
po_ppmx1_aux <- postquant_dm(y = Y, output = out_ppmx1_aux, data = mydata, plot = F)

mytime_ppmx1_dd <- system.time(
  out_ppmx1_dd <- my_dm_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1, similarity = 1, consim=1, calibration=1,
                              similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), modelpriors, mhtune=c(0.5, 0.5),
                              iter=iterations,burn=burnin,thin=thinning))
po_ppmx1_dd <- postquant_dm(y = Y, output = out_ppmx1_dd, data = mydata, plot = F)

mytime_ppmx2_aux <- system.time(
  out_ppmx2_aux <- my_dm_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1, similarity = 1, consim=1, calibration=2,
                    similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), modelpriors, mhtune=c(0.5, 0.5),
                    iter=iterations,burn=burnin,thin=thinning))
po_ppmx2_aux <- postquant_dm(y = Y, output = out_ppmx2_aux, data = mydata, plot = F)

mytime_ppmx2_dd <- system.time(
  out_ppmx2_dd <- my_dm_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1, similarity = 2, consim=1, calibration=2,
                              similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0), modelpriors, mhtune=c(0.5, 0.5),
                              iter=iterations,burn=burnin,thin=thinning))
po_ppmx2_dd <- postquant_dm(y = Y, output = out_ppmx2_dd, data = mydata, plot = F)

tab <- rbind(po_ppm, po_ppmx0_aux, po_ppmx0_dd, po_ppmx1_aux, po_ppmx1_dd,
             po_ppmx2_aux, po_ppmx2_dd)
times <- c(mytime_ppm[3], mytime_ppmx0_aux[3], mytime_ppmx0_dd[3], mytime_ppmx1_aux[3], mytime_ppmx1_dd[3],
           mytime_ppmx2_aux[3], mytime_ppmx2_dd[3])
tab <- cbind(tab[,-4], times)
tab

#postquant_dm(y = Y, output = out_ppm, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx0_aux, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx0_dd, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx1_aux, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx1_dd, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx2_aux, data = mydata, plot = T)
#postquant_dm(y = Y, output = out_ppmx2_dd, data = mydata, plot = T)

