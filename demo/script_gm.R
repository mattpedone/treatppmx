rm(list=ls())

devtools::load_all()

###### STUDIO 2
### Scenario a

KK <- 1#0
res_1a <- matrix(0, KK, 8)
res_2a <- matrix(0, KK, 8)
res_3a <- matrix(0, KK, 8)
res_4a <- matrix(0, KK, 8)
res_5a <- matrix(0, KK, 8)
#set.seed(121)
K=2

modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y))
modelpriors$hP0_L0 <- diag(10, ncol(Y))
modelpriors$hP0_nu0 <- nrow(Y) + 2
modelpriors$hP0_V0 <- diag(10, ncol(Y))

iterations <- 10000
burnin <- 2000
thinning <- 10

nout <- (iterations-burnin)/thinning

##################
## TABELLA 1 DEL SECONDO STUDIO DI SIMULAZIONE
## n = 100, P = 100, Q = 2, K = 2, CC = 5, REUSE = TRUE
## SIMILARITY = AUX, CALIBRATION = COAR, MODEL = 1 (NN)
##################
myppmx <- gcd(n=100, concov = 100, K, alpha = 1, similarity = 1, simparm = 1)
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

for(k in 1:KK){
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 1, consim=1, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_1a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 1, consim=2, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_2a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 2, consim=1, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_3a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 2, consim=2, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_4a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 0,
                       similarity = 2, consim=1, calibration=1,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_5a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)
}
tab1 <- rbind(apply(res_1a, 2, mean), sqrt(apply(res_1a, 2, var)), apply(res_2a, 2, mean),
               sqrt(apply(res_2a, 2, var)), apply(res_3a, 2, mean), sqrt(apply(res_3a, 2, var)),
               apply(res_4a, 2, mean), sqrt(apply(res_4a, 2, var)), apply(res_5a, 2, mean),
               sqrt(apply(res_5a, 2, var)))
colnames(tab1) <- heading
rownames(tab1) <- righe
#save(tab1, file = "tab1.RData")

##################
## TABELLA 2 DEL SECONDO STUDIO DI SIMULAZIONE
## n = 100, P = 100, Q = 2, K = 2, CC = 5, REUSE = TRUE
## SIMILARITY = AUX, CALIBRATION = COAR, MODEL = 2 (NNIG)
##################
myppmx <- gcd(n=100, concov = 100, K, alpha = 1, similarity = 1, simparm = 2)
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

for(k in 1:KK){
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 1, consim=1, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_1a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 1, consim=2, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_2a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 2, consim=1, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_3a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 2, consim=2, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_4a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 0,
                       similarity = 2, consim=1, calibration=1,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_5a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)
}
tab2 <- rbind(apply(res_1a, 2, mean), sqrt(apply(res_1a, 2, var)), apply(res_2a, 2, mean),
              sqrt(apply(res_2a, 2, var)), apply(res_3a, 2, mean), sqrt(apply(res_3a, 2, var)),
              apply(res_4a, 2, mean), sqrt(apply(res_4a, 2, var)), apply(res_5a, 2, mean),
              sqrt(apply(res_5a, 2, var)))
colnames(tab2) <- heading
rownames(tab2) <- righe
#save(tab2, file = "tab2.RData")

##################
## TABELLA 3 DEL SECONDO STUDIO DI SIMULAZIONE
## n = 100, P = 100, Q = 2, K = 2, CC = 5, REUSE = TRUE
## SIMILARITY = DD, CALIBRATION = COAR, MODEL = 1 (NN)
##################
myppmx <- gcd(n=100, concov = 100, K, alpha = 1, similarity = 2, simparm = 1)
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

for(k in 1:KK){
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 1, consim=1, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_1a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 1, consim=2, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_2a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 2, consim=1, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_3a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 2, consim=2, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_4a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 0,
                       similarity = 2, consim=1, calibration=1,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_5a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)
}
tab3 <- rbind(apply(res_1a, 2, mean), sqrt(apply(res_1a, 2, var)), apply(res_2a, 2, mean),
              sqrt(apply(res_2a, 2, var)), apply(res_3a, 2, mean), sqrt(apply(res_3a, 2, var)),
              apply(res_4a, 2, mean), sqrt(apply(res_4a, 2, var)), apply(res_5a, 2, mean),
              sqrt(apply(res_5a, 2, var)))
colnames(tab3) <- heading
rownames(tab3) <- righe
#save(tab3, file = "tab3.RData")

##################
## TABELLA 4 DEL SECONDO STUDIO DI SIMULAZIONE
## n = 100, P = 100, Q = 2, K = 2, CC = 5, REUSE = TRUE
## SIMILARITY = DD, CALIBRATION = COAR, MODEL = 2 (NNIG)
##################
myppmx <- gcd(n=100, concov = 100, K, alpha = 1, similarity = 2, simparm = 2)
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

for(k in 1:KK){
  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 1, consim=1, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_1a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 1, consim=2, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_2a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 2, consim=1, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_3a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                       similarity = 2, consim=2, calibration=2,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_4a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)

  mytime <- system.time(
    out <- my_mvn_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 0,
                       similarity = 2, consim=1, calibration=1,
                       #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                       similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                       modelpriors, mhtune=c(0.5, 0.5),
                       iter=iterations,burn=burnin,thin=thinning))

  res_5a[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)
}
tab4 <- rbind(apply(res_1a, 2, mean), sqrt(apply(res_1a, 2, var)), apply(res_2a, 2, mean),
              sqrt(apply(res_2a, 2, var)), apply(res_3a, 2, mean), sqrt(apply(res_3a, 2, var)),
              apply(res_4a, 2, mean), sqrt(apply(res_4a, 2, var)), apply(res_5a, 2, mean),
              sqrt(apply(res_5a, 2, var)))
colnames(tab4) <- heading
rownames(tab4) <- righe
#save(tab4, file = "tab4.RData")
