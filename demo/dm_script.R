rm(list=ls())

devtools::load_all()

KK <- 1#numero di repliche
res <- matrix(0, KK, 7)
par(mfrow=c(1,1))
K=3#dimensioni
#mydata <- gcd_dm(n_obs = 100, concov = 100, K, similarity = 1, simparm = 1,
#              alpha = 1, m0 = 0, s20 = 1, v = 2, k0 = 10, v0 = 1, plot = F)
#mydata$possmean
#cat("nclus: ", mydata$nclus, "\n")
mydata <- gcd_dm_simpl(100)
Y <- mydata$y

apply(Y, 2, sum)
apply(Y[1:50,], 2, sum)
apply(Y[51:100,], 2, sum)

colors <- c("#ebb678", "#1979a9", "#e07b39", "#69bdd2", "#80391e", "#cce7e8",
            "#1c100b", "#042f66", "#44bcd8")
colors <- colors[mydata$label]
#if(K==2){plot(myppmx$intercept, pch = 16, col = colors)}
#if(K==3){scatterplot3d::scatterplot3d(myppmx$intercept, pch = 16, color = colors, grid=TRUE, box=FALSE)}

heading <- c("nclu_post", "MSE", "lplm", "ARI", "ESS", "nclu_real", "nout")
#righe <- c("Reuse", " ", "No ", " ", "DD Cal.", " ","DD Coa", " ", "PPM", " ")
X <- mydata$X
#plot(X[,1], X[, 2])
modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y))
modelpriors$hP0_L0 <- diag(1, ncol(Y))
modelpriors$hP0_nu0 <- nrow(Y) + 20
modelpriors$hP0_V0 <- diag(100, ncol(Y))

iterations <- 5000
burnin <- 0
thinning <- 1

nout <- (iterations-burnin)/thinning

mytime <- system.time(
  out <- my_dm_ppmx(y = Y, X = X, alpha=.5, CC = 5, reuse = 1, PPMx = 0,
                     similarity = 1, consim=2, calibration=1,
                     #similparam=c(0.0, 10.0, 0.5, 1.0, 10.0, 0.1, 1.0),
                     similparam = c(0.0, 1.0, 0.1, 10.0, 2.0, 0.1, 1.0),
                     modelpriors, mhtune=c(0.5, 0.5),
                     iter=iterations,burn=burnin,thin=thinning))
#res[k,] <- c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout)
#cat("res", c(unlist(postquant(y = Y, output = out, data = myppmx, lab = F, plot = F)), myppmx$nclus, nout), "\n")

postquant_dm(y = Y, output = out, data = mydata, plot = T)#, minbinder = F)
c(mydata$label)
mydata$possmean
apply(out$eta, c(1, 2), mean)
