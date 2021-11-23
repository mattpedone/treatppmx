#### ---- LOOCV ---- ####
rm(list=ls())
set.seed(121, kind = "L'Ecuyer-CMRG")
devtools::load_all()
load("~/Dropbox/PHD/study-treatppmx/data/LGGdata.rda")
#library(treatppmx)
library(parallel)
library(doParallel)
library(mcclust)
library(mcclust.ext)

#name <- c("a1s01.RData")
trtsgn <- c(matchRTComp[,10])
trtsgn[80:158] <- trtsgn[80:158] + 2
npat <- length(trtsgn)

predAPT_all <- matrix(0, nrow = npat, ncol = 9)
nclust_all <- rep(0, 6)
gof_all <- rep(0, 2)
#myres0 <- sellines_all <- vector(mode = "list", length = K)

wk <- c(0, 40, 100)

cor_all <- parallel::detectCores()-1#cores to be allocated
registerDoParallel(cores = cor_all)

X <- data.frame(scale(matchRTComp[,16:38]))
Z <- data.frame(scale(matchRTComp[,c(11,13)]))#data.frame(orgx)#
#aggiusta!
Y <- matrix(0, nrow = nrow(X), ncol = max(as.numeric(matchRTComp[,9])))
for(i in 1:nrow(Y)){
  Y[i, as.numeric(matchRTComp[i,9])] <- 1
}

modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y)); modelpriors$hP0_L0 <- diag(100, ncol(Y))
modelpriors$hP0_nu0 <- ncol(Y) + 20; modelpriors$hP0_V0 <- diag(100, ncol(Y))

n_aux <- 5 # auxiliary variable for Neal's Algorithm 8
vec_par <- c(0.0, 10.0, .5, 1.0, 2.0, 2.0, 0.1)
#double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
iterations <- 5000#0
burnin <- 1000#0
thinning <- 1#0

nout <- (iterations-burnin)/thinning
predAPT <- c()
#sub <- 15
myres <- foreach(sub = 1:npat, .combine = rbind) %dopar%
  {
    out_ppmx <- tryCatch(expr = ppmxct(y = data.matrix(Y[-sub,]), X = data.frame(X[-sub,]),
                                       Xpred = data.frame(X[sub,]), Z = data.frame(Z[-sub,]),
                                       Zpred = data.frame(Z[sub,]), asstreat = trtsgn[-sub], #treatment,
                                       PPMx = 1, cohesion = 2, alpha = 1, sigma = 0.25,
                                       similarity = 2, consim = 1, similparam = vec_par,
                                       calibration = 2, coardegree = 2, modelpriors,
                                       update_hierarchy = T,
                                       hsp = T, iter = iterations, burn = burnin, thin = thinning,
                                       mhtunepar = c(0.05, 0.05), CC = n_aux, reuse = 1, nclu_init = 5), error = function(e){FALSE})

    #number of a cluster, mean, binder &varinf ----
    if(!is.logical(out_ppmx)){
      mc <- apply(out_ppmx$nclu, 1, mean)
      trt <- trtsgn[-sub]
      num_treat <- table(trt)

      cls1 <- t(as.matrix(out_ppmx$label[[1]]))[,c(1:num_treat[1])]
      psm1 <- comp.psm(cls1)
      mc_b1 <- minbinder.ext(psm1);
      mc_vi1 <- minVI(psm1);

      cls2 <- t(as.matrix(out_ppmx$label[[2]]))[,c(1:num_treat[2])]
      psm2 <- comp.psm(cls2)
      mc_b2 <- minbinder.ext(psm2);
      mc_vi2 <- minVI(psm2);

      mc_b <- c(max(mc_b1$cl), max(mc_b2$cl))
      mc_vi <- c(max(mc_vi1$cl), max(mc_vi2$cl))

      #if(sum(apply(out_ppmx$ypred, c(1,2,3), mean)<0)!=0){
      #  #save(out_ppmx$pipred, file=paste0("output/check_", pipred, ".RData"))
      #  #cat("pipred: ", out_ppmx$pipred, "\n")
      #  save(out_ppmx$ypred, file=paste0("output/check_", ypred, ".RData"))
      #}
      #posterior predictive probabilities ----
      A0 <- c(apply(out_ppmx$ypred, c(1,2,3), mean), mc, mc_b, mc_vi, out_ppmx$WAIC, out_ppmx$lpml)}#A0
    ifelse(is.logical(out_ppmx), return(rep(0, 14)), return(A0))
  }

##treatment prediction with utility function
#cat("errori: ", which(rowSums(myres) == 0), "\n")
#myres0[[k]] <- myres
#myres <- myres[complete.cases(myres),]
sellines <- 1:npat#as.vector(which(rowSums(myres) != 0))
#sellines <- 1:npat
A1 <- myres[sellines, 1:3]%*%wk
A2 <- myres[sellines, 4:6]%*%wk
predAPT_all[sellines, 1] <- A1
predAPT_all[sellines, 2] <- A2
myt <- as.numeric(A1 < A2) + 1
predAPT_all[sellines, 3] <- myt
predAPT_all[sellines, 4:9] <- myres[sellines, 1:6]

#NPC
PPMXCUT <- c()
temp <- matrix(0,  nrow = length(sellines), ncol = 6)
temp <- predAPT_all[sellines, 4:9]
npc2 <- function (output, trtsgn, myoutot){
  #K <- dim(output)[3]
  n <- dim(output)[1]
  myctut <- matrix(0, nrow = 3, ncol = 3)
  myctutSum <- NULL
  mycurdata <- output
  myy <- c()
  mypre <- NULL#matrix(0, nrow = nrow(myoutot), ncol = ncol(myoutot))
  pretrt1 <- apply(mycurdata[, 1:3], 1, which.max)
  pretrt2 <- apply(mycurdata[, 4:6], 1, which.max)
  mypreTall <- cbind(pretrt1, pretrt2)
  for (j in 1:n) {
    mypre[j] <- mypreTall[j, trtsgn[j]]
    myy[j] <- match(1, myoutot[j,])
  }
  sts <- table(mypre, myy)
  mysdls <- as.numeric(rownames(sts))
  str1 <- matrix(0, nrow = 3, ncol = 3)
  str1[mysdls, ] <- sts
  myctut <- str1 * diag(3)
  myctutSum <- sum(str1 * diag(3))
  #res <- cbind(myctutSum)
  return(myctutSum)
}
NPC <- npc2(output = temp, trtsgn = trtsgn[sellines], myoutot = Y[sellines,])
NPC

#ESM
myy <- c()
my <- Y[sellines,]
for (j in 1:nrow(my)) {
  myy[j] <- match(1, my[j,]) - 1
}

#ho definito come respondent anche i partial responent
#mytab <- cbind(myass = predAPT_all[sellines,3], rndass = trtsgn[sellines], resp = as.numeric(Y[sellines,1]!=1))
mytab <- cbind(myass = predAPT_all[sellines,3], rndass = trtsgn[sellines], resp = as.numeric(myy>=2))
pred1 <- subset(mytab, mytab[,1]==1)
table1 <- table(pred1[,3],pred1[,2])
pred2 <- subset(mytab, mytab[,1]==2)
table2 <- table(pred2[,3], pred2[,2])
p1 <- sum(table1)/(sum(table1)+sum(table2))
p2 <- sum(table2)/(sum(table1)+sum(table2))

if(length(table1) == 4){
  crt1 <- table1[2,1]/sum(table1[,1])
}
if(length(table1) < 4){
  crt1 <- as.numeric(row.names(table1))
}
if(length(table2) == 4){
  crt2 <- table2[2,2]/sum(table2[,2])
}
if(length(table2) < 4){
  crt2 <- as.numeric(row.names(table2))
}

table1; table2
### summary meaures
esm <- crt1*p1 + crt2*p2 - sum(as.numeric(as.numeric(myy>=2)))/length(myy)
esm

#NCLU
NC <- rbind(mean = apply(myres[sellines, 7:8], 2, mean), sd = apply(myres[sellines, 7:8], 2, sd))
BI <- rbind(mean = apply(myres[sellines, 9:10], 2, mean), sd = apply(myres[sellines, 9:10], 2, sd))
VI <- rbind(mean = apply(myres[sellines, 11:12], 2, mean), sd = apply(myres[sellines, 11:12], 2, sd))

#FIT
gof_all <- rbind(waic = c(mean(myres[sellines, 13]), sd(myres[sellines, 13])),
                 lpml = c(mean(myres[sellines, 14]), sd(myres[sellines, 14])))

#results
cluPPMX <- rbind(NC, BI, VI)
colnames(cluPPMX) <- c("trt 1", "trt 2")
cluPPMX

## Calculate CPO counts
myresults<-cbind(predAPT_all[sellines, 4:9], trtsgn[sellines], as.numeric(paste(matchRTComp[,9])));
ResTable1<-myresults[80:158,c(1:3,8)];              ## received trt1
ResTable2<-myresults[1:79,c(4:6,8)];                ## received trt2

ResTable3<-rbind(ResTable1,ResTable2);

### predicted outcomes;
MoReult<-as.numeric(apply(ResTable3[,1:3],1,which.max))-1;
Cresult<-cbind(ResTable3[,4], MoReult);

table(Cresult[,2],Cresult[,1]);
###  CPO Results
sum(diag(table(Cresult[,2],Cresult[,1])))

#save(resPPMX, file=paste0("output/tuning_scenario3/res_", name))
#save(cluPPMX, file=paste0("output/tuning_scenario3/clu_", name))

#### ---- Single Fold ---- ####
rm(list=ls())
set.seed(121, kind = "L'Ecuyer-CMRG")
devtools::load_all()
load("~/Dropbox/PHD/study-treatppmx/data/LGGdata.rda")
#library(treatppmx)
library(parallel)
library(doParallel)
library(mcclust)
library(mcclust.ext)

#name <- c("a1s01.RData")
trtsgn <- c(matchRTComp[,10])
trtsgn[80:158] <- trtsgn[80:158] + 2
npat <- length(trtsgn)

predAPT_all <- matrix(0, nrow = npat, ncol = 9)
nclust_all <- rep(0, 6)
gof_all <- rep(0, 2)
#myres0 <- sellines_all <- vector(mode = "list", length = K)

wk <- c(0, 40, 100)

cor_all <- parallel::detectCores()-1#cores to be allocated
registerDoParallel(cores = cor_all)

X <- data.frame(scale(matchRTComp[,16:38]))
Z <- data.frame(scale(matchRTComp[,c(11,13)]))#data.frame(orgx)#
#aggiusta!
Y <- matrix(0, nrow = nrow(X), ncol = max(as.numeric(matchRTComp[,9])))
for(i in 1:nrow(Y)){
  Y[i, as.numeric(matchRTComp[i,9])] <- 1
}

modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y)); modelpriors$hP0_L0 <- diag(100, ncol(Y))
modelpriors$hP0_nu0 <- ncol(Y) + 2; modelpriors$hP0_V0 <- diag(1, ncol(Y))

n_aux <- 5 # auxiliary variable for Neal's Algorithm 8
vec_par <- c(0.0, 10.0, .5, 1.0, 2.0, 2.0, 0.1)
#double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
iterations <- 5000
burnin <- 1000
thinning <- 10

nout <- (iterations-burnin)/thinning
sub <- sample(1:npat, 1)
out_ppmx <- ppmxct(y = data.matrix(Y[-sub,]), X = data.frame(X[-sub,]),
                   Xpred = data.frame(X[sub,]), Z = data.frame(Z[-sub,]),
                   Zpred = data.frame(Z[sub,]), asstreat = trtsgn[-sub], #treatment,
                   PPMx = 1, cohesion = 2, alpha = 1, sigma = 0.25,
                   similarity = 2, consim = 1, similparam = vec_par,
                   calibration = 2, coardegree = 2, modelpriors,
                   update_hierarchy = T,
                   hsp = T, iter = iterations, burn = burnin, thin = thinning,
                   mhtunepar = c(0.05, 0.05), CC = n_aux, reuse = 1, nclu_init = 5)

#number of a cluster, mean, binder &varinf
mc <- apply(out_ppmx$nclu, 1, mean)
trt <- trtsgn[-sub]
num_treat <- table(trt)

cls1 <- t(as.matrix(out_ppmx$label[[1]]))[,c(1:num_treat[1])]
psm1 <- comp.psm(cls1)
mc_b1 <- minbinder.ext(psm1);
mc_vi1 <- minVI(psm1);

cls2 <- t(as.matrix(out_ppmx$label[[2]]))[,c(1:num_treat[2])]
psm2 <- comp.psm(cls2)
mc_b2 <- minbinder.ext(psm2);
mc_vi2 <- minVI(psm2);

mc_b <- c(max(mc_b1$cl), max(mc_b2$cl)); mc_b
mc_vi <- c(max(mc_vi1$cl), max(mc_vi2$cl)); mc_vi

# Traceplot for the number of clusters
df <- data.frame(t(out_ppmx$nclu))
colnames(df) <- c("t1", "t2")
df <- cbind(Index = as.numeric(row.names(df)), df)
df <- reshape2::melt(df, id.vars="Index")
ggplot2::ggplot(df, aes(x = Index, y = value, col = variable)) + geom_line() + theme_classic()

## Eta posterior
#moda <- function(v) {
#  tmp <- unique(v)
#  v[which.max(tabulate(match(v, tmp)))]
#}
#
#df <- data.frame(t(out_ppmx$nclu))
#mode <- apply(df, 2, moda)
#
#idx1 <- df[,1] == mode[1]
#idx2 <- df[,2] == mode[2]
#
#nmax <- length(unlist(out_ppmx$eta[1,1]))/3
#
#eta1 <- array(unlist(out_ppmx$eta[idx1, 1]), dim = c(3, nmax, sum(idx1)))[,1:mode[1],]
#
#plot(density(eta1[1,1,]))
#plot(density(eta1[2,1,]))
#plot(density(eta1[3,1,]))
#
#plot(density(eta1[1,2,]))
#plot(density(eta1[2,2,]))
#plot(density(eta1[3,2,]))
#
#plot(density(eta1[1,3,]))
#plot(density(eta1[2,3,]))
#plot(density(eta1[3,3,]))
#
#plot(density(eta1[1,4,]))
#plot(density(eta1[2,4,]))
#plot(density(eta1[3,4,]))


#posterior predictive probabilities
A0 <- c(apply(out_ppmx$ypred, c(1,2,3), mean), mc, mc_b, mc_vi, out_ppmx$WAIC, out_ppmx$lpml)
A0

#### ---- Check RNG ---- ####

out <- treatppmx::ran_mvnorm(rep(0, 3), c(diag(1, 3)), 3)
for(i in 1:10000){
  out <- cbind(out, treatppmx::ran_mvnorm(rep(0, 3), c(diag(1, 3)), 3))
}
apply(out, 1, mean)
apply(out, 1, var)

out <- rmvnorm(1, rep(0, 3), sigma = diag(1, 3))
for(i in 1:10000){
  out <- rbind(out, rmvnorm(1, rep(0, 3), sigma = diag(1, 3)))
}
apply(out, 2, mean)
apply(out, 2, var)



prob <- runif(200)
prob <- prob/sum(prob) # standardize the probabilities
size <- 500
n <- 20

set.seed(10)
sim_r <- c(rmvnorm(1, rep(0, 3), sigma = diag(10, 3)))
set.seed(10)
sim_rcpp <- treatppmx::ran_mvnorm(rep(0, 3), c(diag(10, 3)), 3)
all.equal(sim_r, sim_rcpp)

microbenchmark::microbenchmark(
  "r" = c(rmvnorm(1, rep(0, 3), sigma = diag(10, 3))),
  "rcpp" = c(treatppmx::ran_mvnorm(rep(0, 3), c(diag(10, 3)), 3))
)


S <- toeplitz((10:1)/10)
set.seed(11)
R <- rWishart(10000, 20, S)
dim(R)  #  10 10  1000
mR <- apply(R, 1:2, mean)  # ~= E[ Wish(S, 20) ] = 20 * S
stopifnot(all.equal(mR, 20*S, tolerance = .01))

L0l <- list()
for(k in 1:10000){
L0v <- treatppmx::ran_iwish(nu = 20, S, dim = 10)
L0 <- matrix(0, 10, 10)

for(i in 0:9){
  for(j in 0:9){
    L0[i+1, j+1] = L0v[(i*10 + j)+1]
  }
}

R[,,k] <- L0}
mR <- apply(R, 1:2, mean)  # ~= E[ Wish(S, 20) ] = 20 * S
stopifnot(all.equal(mR, 20*S, tolerance = .01))
