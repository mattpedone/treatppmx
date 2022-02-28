devtools::load_all()

k=1 #sample(1:30, 1)
load("~/Dropbox/PHD/study-treatppmx/data/scenario3.rda")
print(k)
X <- data.frame(mydata)
Z <- data.frame(cbind(myz2, myz3))
Y <- mytot[,,k]
idx <- sort(sample(1:nrow(Y), 1, replace = F))
print(idx)
wk <- c(0, 40, 100)
df <- data.frame(myprob[[2]]%*%(wk)-myprob[[1]]%*%(wk))
colnames(df) <- c("Utility")
df <- cbind(Iteration = as.numeric(row.names(df)), df)
ggplot2::ggplot(df, aes(x = Iteration, y = Utility)) + geom_line() + theme_classic()

treattest <- trtsgn[idx]
trt <- trtsgn[-idx]
Xtest <- data.frame(X[idx,])
X <- data.frame(X[-idx,])
Ytest <- data.matrix(Y[idx,])
Y <- data.matrix(Y[-idx,])
Ztest <- data.frame(Z[idx,])
Z <- data.frame(Z[-idx,])
nobs <- nrow(Y)

optrt <- as.numeric(myprob[[2]][idx,]%*%wk > myprob[[1]][idx,]%*%wk)+1; optrt

modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y)); modelpriors$hP0_L0 <- diag(10, ncol(Y))
modelpriors$hP0_nu0 <- ncol(Y) + 2; modelpriors$hP0_V0 <- diag(1.0, ncol(Y))

n_aux <- 5
vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
#double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
iterations <- 2000
burnin <- 0
thinning <- 1

nout <- (iterations-burnin)/thinning
time_ppmx <- system.time(
  out_ppmx <- ppmxct(y = Y, X = X, Xpred = Xtest, Z = Z, Zpred = Ztest,
                     asstreat = trt, PPMx = 1, kappa = c(1, 5, 10, 1), sigma = c(0.005, .5, 20),
                     CC = n_aux, cohesion = 2, similarity = 2, consim = 2,
                     calibration = 2, coardegree = 2, similparam = vec_par,
                     modelpriors = modelpriors, iter = iterations,
                     burn = burnin, thin = thinning, nclu_init = 10))
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

# Posterior frequency for (\kappa, \sigma) ----

par(mfrow=c(2,1))
hist(out_ppmx$sigmangg[1,], breaks = 10)
hist(out_ppmx$sigmangg[2,], breaks = 10)
plot(out_ppmx$sigmangg[1,], type ="l")
plot(out_ppmx$sigmangg[2,], type ="l")
table(out_ppmx$sigmangg[1,])
table(out_ppmx$sigmangg[2,])

hist(out_ppmx$kappangg[1,], breaks = 10)
hist(out_ppmx$kappangg[2,], breaks = 10)
plot(out_ppmx$kappangg[1,], type ="l")
plot(out_ppmx$kappangg[2,], type ="l")
table(out_ppmx$kappangg[1,])
table(out_ppmx$kappangg[2,])

P <- table(out_ppmx$sigmangg[1,], out_ppmx$kappangg[1,])
Pm <- reshape::melt(P)
ggplot2::ggplot(Pm, aes(Var.1, Var.2, fill=value)) + geom_tile() +
  ggplot2::geom_text(aes(label=value),colour="white")

P <- table(out_ppmx$sigmangg[2,], out_ppmx$kappangg[2,])
Pm <- reshape::melt(P)
ggplot2::ggplot(Pm, aes(Var.1, Var.2, fill=value)) + geom_tile() +
  ggplot2::geom_text(aes(label=value),colour="white")

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
A0 <- apply(out_ppmx$ypred, c(1,2,3), mean, na.rm=TRUE);#A0
A0 <- c(apply(out_ppmx$pipred, c(1,2,3), median, na.rm=TRUE))#, mc, mc_b, mc_vi, out_ppmx$WAIC, out_ppmx$lpml)

#treatmente prediction with utility function ----
optrt <- as.numeric(myprob[[2]][idx,]%*%wk > myprob[[1]][idx,]%*%wk)+1
#predtrt <- as.numeric(A0[,,2]%*%wk > A0[,,1]%*%wk)+1
predtrt <- as.numeric(A0[4:6]%*%wk > A0[1:3]%*%wk)+1
print(optrt); print(predtrt)

sum(optrt==predtrt)/length(predtrt)

#posterior distribution of predictive utility
ns <- dim(out_ppmx$pipred)[4]

dpu <- matrix(0, ns, 2)
for(i in 1:ns){
  dpu[i,] <- apply(out_ppmx$pipred[,,,i]*wk, 2, sum)
}

par(mfrow=c(2, 1))
mymean <- apply(out_ppmx$pipred, c(1,2,3), mean, na.rm=TRUE); sum(mymean[,,1] * wk) - sum(mymean[,,2] * wk)
mymedian <- apply(out_ppmx$pipred, c(1,2,3), median, na.rm=TRUE); sum(mymedian[,,1] * wk) - sum(mymedian[,,2] * wk)
#plot(density(dpu[,1]), ylim = c(0, .055))
hist(dpu[,1], breaks = 20)
abline(v = mymean[1:3]%*%wk, col = "red")
abline(v = mymedian[1:3]%*%wk, col = "blue")
#plot(density(dpu[,2]), ylim = c(0, .055))
hist(dpu[,2], breaks = 20)
abline(v = mymean[4:6]%*%wk, col = "red")
abline(v = mymedian[4:6]%*%wk, col = "blue")

