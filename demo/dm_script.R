rm(list=ls())

devtools::load_all()

###### STUDIO 1
### Scenario b
#set.seed(121)
KK <- 1#numero di repliche
res <- matrix(0, KK, 7)
par(mfrow=c(1,1))
K=3#dimensioni
myppmx <- gcd_dm(n_obs = 100, concov = 20, K, similarity = 1, simparm = 1,
              alpha = 1, m0 = 0, s20 = 1, v = 2, k0 = 10, v0 = 1, plot = F)
cat("nclus: ", myppmx$nclus, "\n")
Y <- myppmx$y
colors <- c("#ebb678", "#1979a9", "#e07b39", "#69bdd2", "#80391e", "#cce7e8",
            "#1c100b", "#042f66", "#44bcd8")
colors <- colors[myppmx$label]
#if(K==2){plot(Y, pch = 16, col = colors)}
#if(K==3){scatterplot3d(Y, pch = 16, color = colors, grid=TRUE, box=FALSE)}

heading <- c("nclu_post", "MSE", "lplm", "ARI", "ESS", "nclu_real", "nout")
#righe <- c("Reuse", " ", "No ", " ", "DD Cal.", " ","DD Coa", " ", "PPM", " ")
X <- myppmx$X

modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y))
modelpriors$hP0_L0 <- diag(10, ncol(Y))
modelpriors$hP0_nu0 <- nrow(Y) + 2
modelpriors$hP0_V0 <- diag(10, ncol(Y))

iterations <- 500
burnin <- 10#2000
thinning <- 1#10

nout <- (iterations-burnin)/thinning

mytime <- system.time(
  out <- my_dm_ppmx(y = Y, X = X, alpha=1, CC = 5, reuse = 1, PPMx = 1,
                     similarity = 1, consim=1, calibration=1,
                     #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                     similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                     modelpriors, mhtune=c(0.5, 0.5),
                     iter=iterations,burn=burnin,thin=thinning))

#res[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)
#cat("res", c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout), "\n")
postquant_dm(y = Y, output = out, data = myppmx, plot = F)
