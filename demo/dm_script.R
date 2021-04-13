rm(list=ls())

#devtools::load_all()
library(mvtnorm)
library(mcclust)
library(mclust)
library(coda)
library(mcclust.ext)
Rcpp::sourceCpp(file = "src/utils.cpp")
Rcpp::sourceCpp(file = "src/dm_ppmx.cpp")
source(file = "R/dm_ppmx.R")
source(file = "R/rppmx.R")

KK <- 1#numero di repliche
res <- matrix(0, KK, 7)
par(mfrow=c(1,1))
K = 4#dimensioni
mydata <- gcd_dm(n_obs = 100, concov = 2, K, similarity = 1, simparm = 1,
              alpha = 1, m0 = 0, s20 = 1, v = 2, k0 = 10, v0 = 1, plot = F)
mydata$possmean
cat("nclus: ", mydata$nclus, "\n")
#load("data/example1.RData")

#mydata <- gcd_dm_simpl()
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
modelpriors$hP0_nu0 <- nrow(Y) + 2
modelpriors$hP0_V0 <- diag(100, ncol(Y))

alpha_DP <- 1
n_aux <- 5
vec_par <- c(0.0, 10.0, .5, 1.0, 2.0, 2.0, 0.1)
#double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
mhtune=c(0.5, 0.5)
iterations <- 10000
burnin <- 0
thinning <- 1

nout <- (iterations-burnin)/thinning

# PPM w/o Reuse
#time_ppm_nr <- system.time(
#  out_ppm_nr <- my_dm_ppmx(y = Y, X = X, alpha=alpha_DP, CC = n_aux, reuse = 0,
#                           PPMx = 0, similarity = 1, consim=1, calibration=0,
#                           similparam = vec_par, modelpriors, update_hierarchy = F,
#                           iter=iterations,burn=burnin,thin=thinning))
#ppm_nr <- postquant_dm(y = Y, output = out_ppm_nr, data = mydata, plot = T)

# PPM w Reuse
time_ppm <- system.time(
  out_ppm <- my_dm_ppmx(y = Y, X = X, alpha=alpha_DP, CC = n_aux, reuse = 1,
                        PPMx = 0, similarity = 1, consim=1, calibration=0,
                        similparam = vec_par, modelpriors, update_hierarchy = F,
                        iter=iterations,burn=burnin,thin=thinning))
ppm <- postquant_dm(y = Y, output = out_ppm, data = mydata, plot = F)

# PPMx No Calibration Auxiliary similarity w/o Reuse
#time_ppmx0_aux_nr <- system.time(
#  out_ppmx0_aux_nr <- my_dm_ppmx(y = Y, X = X, alpha=alpha_DP, CC = n_aux, reuse = 0,
#                              PPMx = 1, similarity = 1, consim=1, calibration=0,
#                              similparam = vec_par, modelpriors, update_hierarchy = F,
#                              iter=iterations,burn=burnin,thin=thinning))
#ppmx0_aux_nr <- postquant_dm(y = Y, output = out_ppmx0_aux_nr, data = mydata, plot = F)

# PPMx No Calibration Auxiliary Similarity w Reuse

time_ppmx0_aux <- system.time(
  out_ppmx0_aux <- my_dm_ppmx(y = Y, X = X, alpha=alpha_DP, CC = n_aux, reuse = 1,
                              PPMx = 1, similarity = 1, consim=1, calibration=0,
                              similparam = vec_par, modelpriors,  update_hierarchy = F,
                              iter=iterations,burn=burnin,thin=thinning))
ppmx0_aux <- postquant_dm(y = Y, output = out_ppmx0_aux, data = mydata, plot = F)

out_ppmx0_aux$acc_rate_eta
apply(out_ppmx0_aux$eta, c(1,2), mean)
mydata$possmean

# PPMx No Calibration Double Dipper similarity w/o Reuse
#time_ppmx0_dd_nr <- system.time(
#  out_ppmx0_dd_nr <- my_dm_ppmx(y = Y, X = X, alpha=alpha_DP, CC = n_aux, reuse = 0,
#                             PPMx = 1, similarity = 2, consim=1, calibration=0,
#                             similparam = vec_par, modelpriors, update_hierarchy = F,
#                             iter=iterations,burn=burnin,thin=thinning))
#ppmx0_dd_nr <- postquant_dm(y = Y, output = out_ppmx0_dd_nr, data = mydata, plot = F)

# PPMx No Calibration Double Dipper similarity w Reuse
time_ppmx0_dd <- system.time(
  out_ppmx0_dd <- my_dm_ppmx(y = Y, X = X, alpha=alpha_DP, CC = n_aux, reuse = 1,
                             PPMx = 1, similarity = 2, consim=1, calibration=0,
                             similparam = vec_par, modelpriors, update_hierarchy = F,
                             iter=iterations,burn=burnin,thin=thinning))
ppmx0_dd <- postquant_dm(y = Y, output = out_ppmx0_dd, data = mydata, plot = F)

# PPMx Calibrated Auxiliary similarity w/o Reuse
#time_ppmx1_aux_nr <- system.time(
#  out_ppmx1_aux_nr <- my_dm_ppmx(y = Y, X = X, alpha=alpha_DP, CC = n_aux, reuse = 0,
#                              PPMx = 1, similarity = 1, consim=1, calibration=1,
#                              similparam = vec_par, modelpriors,
#                              iter=iterations,burn=burnin,thin=thinning))
#ppmx1_aux_nr <- postquant_dm(y = Y, output = out_ppmx1_aux_nr, data = mydata, plot = F)

# PPMx Calibrated Auxiliary similarity w Reuse
time_ppmx1_aux <- system.time(
  out_ppmx1_aux <- my_dm_ppmx(y = Y, X = X, alpha=alpha_DP, CC = n_aux, reuse = 1,
                              PPMx = 1, similarity = 1, consim=1, calibration=1,
                              similparam = vec_par, modelpriors, update_hierarchy = F,
                              iter=iterations,burn=burnin,thin=thinning))
ppmx1_aux <- postquant_dm(y = Y, output = out_ppmx1_aux, data = mydata, plot = F)

# PPMx Calibrated Double Dipper similarity w/o Reuse
#time_ppmx1_dd_nr <- system.time(
#  out_ppmx1_dd_nr <- my_dm_ppmx(y = Y, X = X, alpha=alpha_DP, CC = n_aux, reuse = 0,
#                             PPMx = 1, similarity = 2, consim=1, calibration=1,
#                             similparam = vec_par, modelpriors,
#                             iter=iterations,burn=burnin,thin=thinning))
#ppmx1_dd_nr <- postquant_dm(y = Y, output = out_ppmx1_dd_nr, data = mydata, plot = F)

# PPMx Calibrated Double Dipper similarity w Reuse
time_ppmx1_dd <- system.time(
  out_ppmx1_dd <- my_dm_ppmx(y = Y, X = X, alpha=alpha_DP, CC = n_aux, reuse = 1,
                             PPMx = 1, similarity = 2, consim=1, calibration=1,
                             similparam = vec_par, modelpriors, update_hierarchy = F,
                             iter=iterations,burn=burnin,thin=thinning))
ppmx1_dd <- postquant_dm(y = Y, output = out_ppmx1_dd, data = mydata, plot = F)

# PPMx Coarsened Auxiliary similarity w/o Reuse
#time_ppmx2_aux_nr <- system.time(
#  out_ppmx2_aux_nr <- my_dm_ppmx(y = Y, X = X, alpha=alpha_DP, CC = n_aux, reuse = 0,
#                                 PPMx = 1, similarity = 1, consim=1, calibration=2,
#                                 similparam = vec_par, modelpriors, update_hierarchy = F,
#                                 iter=iterations,burn=burnin,thin=thinning))
#ppmx2_aux_nr <- postquant_dm(y = Y, output = out_ppmx2_aux_nr, data = mydata, plot = F)


# PPMx Coarsened Auxiliary similarity w Reuse
time_ppmx2_aux <- system.time(
  out_ppmx2_aux <- my_dm_ppmx(y = Y, X = X, alpha=alpha_DP, CC = n_aux, reuse = 1,
                              PPMx = 1, similarity = 1, consim=1, calibration=2,
                              similparam = vec_par, modelpriors, update_hierarchy = F,
                              iter=iterations,burn=burnin,thin=thinning))
ppmx2_aux <- postquant_dm(y = Y, output = out_ppmx2_aux, data = mydata, plot = F)

# PPMx Coarsened Double Dipper similarity w/o Reuse
#time_ppmx2_dd_nr <- system.time(
#  out_ppmx2_dd_nr <- my_dm_ppmx(y = Y, X = X, alpha=alpha_DP, CC = n_aux, reuse = 0,
#                             PPMx = 1, similarity = 2, consim=1, calibration=2,
#                             similparam = vec_par, modelpriors, update_hierarchy = F,
#                             iter=iterations,burn=burnin,thin=thinning))
#ppmx2_dd_nr <- postquant_dm(y = Y, output = out_ppmx2_dd_nr, data = mydata, plot = F)

# PPMx Coarsened Double Dipper similarity w Reuse
time_ppmx2_dd <- system.time(
  out_ppmx2_dd <- my_dm_ppmx(y = Y, X = X, alpha=alpha_DP, CC = n_aux, reuse = 1,
                             PPMx = 1, similarity = 2, consim=1, calibration=2,
                             similparam = vec_par, modelpriors, update_hierarchy = F,
                             iter=iterations,burn=burnin,thin=thinning))
ppmx2_dd <- postquant_dm(y = Y, output = out_ppmx2_dd, data = mydata, plot = F)

#tab <- rbind(ppm_nr, ppm, ppmx0_aux_nr, ppmx0_aux, ppmx0_dd_nr, ppmx0_dd,
#             ppmx1_aux_nr, ppmx1_aux, ppmx1_dd_nr, ppmx1_dd, ppmx2_aux_nr,
#             ppmx2_aux, ppmx2_dd_nr, ppmx2_dd)
#time <- c(time_ppm_nr[3], time_ppm[3], time_ppmx0_aux_nr[3], time_ppmx0_aux[3],
#          time_ppmx0_dd_nr[3], time_ppmx0_dd[3], time_ppmx1_aux_nr[3], time_ppmx1_aux[3],
#          time_ppmx1_dd_nr[3], time_ppmx1_dd[3], time_ppmx2_aux_nr[3],
#          time_ppmx2_aux[3], time_ppmx2_dd_nr[3], time_ppmx2_dd[3])

tab <- rbind(ppm, ppmx0_aux, ppmx0_dd, ppmx1_aux, ppmx1_dd, ppmx2_aux, ppmx2_dd)
time <- c(time_ppm[3], time_ppmx0_aux[3], time_ppmx0_dd[3], time_ppmx1_aux[3],
          time_ppmx1_dd[3], time_ppmx2_aux[3], time_ppmx2_dd[3])
tab <- cbind(tab[,-4], time)
tab

sum(unlist(tab[,4]))/60

postquant_dm(y = Y, output = out_ppm, data = mydata, plot = T)
postquant_dm(y = Y, output = out_ppmx0_aux, data = mydata, plot = T)
postquant_dm(y = Y, output = out_ppmx0_dd, data = mydata, plot = T)
postquant_dm(y = Y, output = out_ppmx1_aux, data = mydata, plot = T)
postquant_dm(y = Y, output = out_ppmx1_dd, data = mydata, plot = T)
postquant_dm(y = Y, output = out_ppmx2_aux, data = mydata, plot = T)
postquant_dm(y = Y, output = out_ppmx2_dd, data = mydata, plot = T)

tab

