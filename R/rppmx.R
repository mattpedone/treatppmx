ran_ppmx <- function(X=NULL, similarity = 1, simparm = 1, alpha=1, m0=0, s20=1,v=2, k0=10, v0=1){

  out <- NULL

  nobs <- nrow(X)
  if(!is.data.frame(X)) X <- data.frame(X)

  classes <- sapply(X, class)
  catvars <- classes %in% c("factor","character")

  # standardize continuous covariates
  if(sum(!catvars) > 0){
    xcon <- apply(X[,!catvars, drop=FALSE], 2, scale)
    ncon <- ncol(xcon)
  }else{
    xcon <- cbind(rep(0,nobs));
    ncon <- 0
  }

  # Function that relabels categorical variables to begin with 0
  relab <- function(x) as.numeric(as.factor(as.character(x))) - 1

  if(sum(catvars) > 0){
    # Change the factors or characters into integers with category starting at 0.
    xcat <- apply(X[, catvars,drop=FALSE], 2, relab)
    Cvec <- apply(xcat,2,function(x)length(unique(x)))
    ncat <- ncol(xcat)
  }else{
    xcat <- cbind(rep(0,nobs));
    Cvec <- 0
    ncat <- 0
  }
  nk <- 1
  nh <- rep(0, nobs)
  Si <- rep(0, nobs)

  dirweights <- rep(0.1, length=max(Cvec));
  #N <- m

  #cat("xcon", as.vector(t(xcon)), "\n")
  #cat("xcat", as.vector(t(xcat)), "\n")
  #cat("Cvec", as.vector(t(Cvec)), "\n")
  #cat("m0", as.double(m0), "\n")
  #cat("k0", as.double(k0), "\n")
  #cat("v0", as.double(v0), "\n")
  #cat("s20", as.double(s20), "\n")
  #cat("v", as.double(v), "\n")
  #cat("dw", as.vector(dirweights), "\n")

  Cout <- ranppmx(as.integer(nobs), as.integer(similarity), as.integer(simparm), as.double(alpha),
                as.integer(ncon), as.integer(ncat), as.vector(t(xcon)),
                as.vector(t(xcat)), as.vector(t(Cvec)), as.double(m0),
                as.double(k0), as.double(v0), as.double(s20),
                as.double(v), as.vector(dirweights))


  label <- Cout$cluster_label
  nclus <- Cout$nclus
  nj <- Cout$nj
  out$label <- label
  out$nclus <- nclus
  #cat("check", nclus, "\n")
  out$nj <- nj[1:nclus]
  return(out)
}

genera_dati <- function(n = 100, P = 2, Q = 2){
  discretize <- function(vec){
    vec[which(vec<0)] <- 0
    vec[which(vec>0)] <- 1
    #vec <- as.factor(vec)
    return(vec)
  }
  ncov <- Q + P
  ##prima ci stanno le binarie
  XX <- matrix(0, n, (ncov))

  dim = 3
  K = dim - 1 #numero di componenti di una mistura -1
  tab <- matrix(0, dim, ncov)
  tab[1, ] <- runif(ncov, -5, 5)

  for(k in 1:K){
    cond = FALSE
    while(cond == FALSE){
      vec <- runif(ncov, -10, 10)
      cond <- all(dist(rbind(tab[c(1:k),], vec))/ncov>2)
    }
    tab[k+1,] <- vec
  }
  #cat("tab: ", "\n", tab, "\n")
  #stick_breaking_process = function(num_weights, alpha) {
  #  betas = rbeta(num_weights, 1, alpha)
  #  remaining_stick_lengths = c(1, cumprod(1 - betas))[1:num_weights]
  #  weights = remaining_stick_lengths * betas
  #  weights
  #}
  #tab[1,] <- stick_breaking_process(ncov, .1)*10
  #tab[2,] <- sort(stick_breaking_process(ncov, .1)*10)
  #tab[3,] <- -tab[1,]


  pro <- c(0.2,0.5,0.3)
  clusterlabel <- c()
  for(i in 1:n){
    u <- runif(1)
    if(u<pro[1]){
      if(Q!=0){XX[i, c(1:Q)] <- discretize(rmvnorm(1, tab[1, c(1:Q)]))}
      if(P!=0){XX[i, c((Q+1):ncov)] <- rmvnorm(1, tab[1, c((Q+1):ncov)])}
    }
    else{
      if(u<(pro[1]+pro[2])){
        if(Q!=0){XX[i, c(1:Q)] <- discretize(rmvnorm(1, tab[2, c(1:Q)]))}
        if(P!=0){XX[i, c((Q+1):ncov)] <- rmvnorm(1, tab[2, c((Q+1):ncov)])}
      }
      else{
        if(Q!=0){XX[i, c(1:Q)] <- discretize(rmvnorm(1, tab[3, c(1:Q)]))}
        if(P!=0){XX[i, c((Q+1):ncov)] <- rmvnorm(1, tab[3, c((Q+1):ncov)])}
      }
    }
  }
  #XX[, seq(1:Q)] <- factor(XX[, seq(1:Q)])
  #XX <- data.frame(XX)
  #X <- XX[, -seq(1:Q)]
  return(list(XX=data.frame(XX)))
}

#' gendata
#'
#' @export
#'
gcd <- function(n_obs, concov = 2, K = 2, similarity = 1, simparm = 1,
                alpha = 1, m0 = 0, s20 = 1, v = 2, k0 = 10, v0 = 1, plot = F){
  myinfopart <- NULL
  #concov = 2
  #n_obs = 100
  n = n_obs
  d <- genera_dati(n=n, P=concov)

  X <- d$XX

  X$X1 <- as.factor(X$X1)
  X$X2 <- as.factor(X$X2)

  myinfopart$X <- X

  myppmx <- ran_ppmx(X=X, similarity, simparm, alpha, m0, s20,v, k0, v0)
  myinfopart$label <- myppmx$label
  myinfopart$nclus <- myppmx$nclus
  myinfopart$nj <- myppmx$nj

  possmean <- matrix(0, myppmx$nclus, K)
  for(i in 1:myppmx$nclus){
    possmean[i,] <- rmvnorm(1, rep(0, K), sigma = diag(50, K), method="svd")
  }
  myinfopart$possmean <- possmean
  Y <- matrix(0, n, K)
  for(i in 1:n){
    Y[i, ] <- rmvnorm(1, possmean[myppmx$label[i],], sigma = diag(.25, K), method="svd")
  }
  myinfopart$y <- Y

  if(plot==T){
    par(mfrow=c(1, 2))
    plot(possmean)
    plot(Y)
  }
  return(myinfopart)
}

#' gendata dm
#'
#' @export
#'
gcd_dm <- function(n_obs, concov = 2, K, similarity = 1, simparm = 1,
                alpha = 1, m0 = 0, s20 = 1, v = 2, k0 = 10, v0 = 1, plot = F){
  myinfopart <- NULL
  #concov = 2
  #n_obs = 100
  n = n_obs
  d <- genera_dati(n=n, P=concov)

  X <- d$XX

  X$X1 <- as.factor(X$X1)
  X$X2 <- as.factor(X$X2)

  myinfopart$X <- X

  myppmx <- ran_ppmx(X=X, similarity, simparm, alpha, m0, s20,v, k0, v0)
  myinfopart$label <- myppmx$label
  myinfopart$nclus <- myppmx$nclus
  myinfopart$nj <- myppmx$nj

  possmean <- matrix(0, myppmx$nclus, K)
  for(i in 1:myppmx$nclus){
    possmean[i,] <- rmvnorm(1, rep(0, K), sigma = diag(10, K), method="svd")
  }
  #cat("possmean", possmean, "\n")
  myinfopart$possmean <- possmean
  Y <- matrix(0, n, K)
  intercept <- matrix(0, n, K)

  for(i in 1:n){
    #intercept[i,] <- rmvnorm(1, possmean[myppmx$label[i], ], sigma = diag(.25, K), method = "svd")
    intercept[i,] <- possmean[myppmx$label[i], ]
  }
  myinfopart$intercept <- intercept

  for(i in 1:n){
    thisrow = as.vector(exp(intercept[i,]))# %*% XX[ii, ]))
    pi = thisrow/sum(thisrow)
    Y[i, ] = rmultinom(1, 1, pi)
  }

  myinfopart$y <- Y

  if(plot==T){
    par(mfrow=c(1, 2))
    plot(possmean)
    plot(Y)
  }
  return(myinfopart)
}

#' Scenario 1
#'
#' @export
#'
scenario1 <- function(){

  dat <- NULL
  n1 <- 85
  n2 <- 30
  n3 <- 85

  n <- n1 + n2 + n3

  x <- matrix(0, nrow = n, ncol = 4)
  label <- c()

  for(i in 1:n){
    if(i <= (n1)){
      x[i,1] <- rnorm(1, -3, sqrt(.5))
      x[i,2] <- rnorm(1, 3, sqrt(.5))
      x[i,3] <- rbinom(1, 1, .1)
      x[i,4] <- rbinom(1, 1, .1)
      label[i] <- 1
    }
    if((n1 < i) & (i <= (n1+n2))){
      x[i,1] <- rnorm(1, 0, sqrt(.5))
      x[i,2] <- rnorm(1, 0, sqrt(.5))
      x[i,3] <- rbinom(1, 1, .1)
      x[i,4] <- rbinom(1, 1, .9)
      label[i] <- 2
    }
    if((i > (n1+n2))){
      x[i,1] <- rnorm(1, 3, sqrt(.5))
      x[i,2] <- rnorm(1, -3, sqrt(.5))
      x[i,3] <- rbinom(1, 1, .1)
      x[i,4] <- rbinom(1, 1, .1)
      label[i] <- 3
    }
  }

  intercept <- x
  Y <- matrix(0, n, 4)

  for(i in 1:n){
    thisrow = as.vector(exp(intercept[i,]))# %*% XX[ii, ]))
    pi = bayess::rdirichlet(n = 1, par = thisrow)
    #pi = thisrow/sum(thisrow)
    Y[i, ] = rmultinom(1, 1, pi)
  }

  x <- data.frame(x)

  x$X3 <- as.factor(x$X3)
  x$X4 <- as.factor(x$X4)

  dat$X <- x
  dat$label <- label
  dat$intercept <- intercept
  dat$y <- Y
  dat$nclus <- 3
  return(dat)
}

#' Scenario 1 in paper Avis modified
#'
#' @export
#'
scenario2 <- function(){

  dat <- NULL
  n1 <- 80
  n2 <- 40#10
  n3 <- 80

  n <- n1 + n2 + n3

  Q = 2
  z <- matrix(0, nrow = n, Q)
  mz <- sample(seq(-2, 2, length.out = 100), Q)
  for(q in 1:Q){
    z[,q] <- rnorm(nrow(z), mz[q], .5)#
  }
  x <- matrix(0, nrow = n, ncol = 4)

  label <- c()

  for(i in 1:n){
    if(i <= (n1)){
      x[i,1] <- rnorm(1, -3, sqrt(.5))
      x[i,2] <- rnorm(1, 3, sqrt(.5))
      x[i,3] <- rbinom(1, 1, .5)
      x[i,4] <- rbinom(1, 1, .1)
      label[i] <- 1
    }
    if((n1 < i) & (i <= (n1+n2))){
      x[i,1] <- rnorm(1, 0, sqrt(.5))
      x[i,2] <- rnorm(1, 0, sqrt(.5))
      x[i,3] <- rbinom(1, 1, .5)
      x[i,4] <- rbinom(1, 1, .9)
      label[i] <- 2
    }
    if((i > (n1+n2))){
      x[i,1] <- rnorm(1, 3, sqrt(.5))
      x[i,2] <- rnorm(1, -3, sqrt(.5))
      x[i,3] <- rbinom(1, 1, .5)
      x[i,4] <- rbinom(1, 1, .1)
      label[i] <- 3
    }
  }
  beta1 <- c(3, 2, 1, 0)/2
  beta2 <- c(-2, -2, 1, 3)/2
  beta3 <- c(3, -2, -1, -1)/2

  theta <- matrix(0, 4, Q)
  #st = 0
  #low_side = .5#beta_min
  #high_side = 1.5#beta_max
  #n_relevant_taxa = 3
  #n_relevant_x = 2
  #if (n_relevant_taxa != 1) {
  #  # warning if the lengths don't match
  #  coef = suppressWarnings(seq(low_side, high_side, len = n_relevant_taxa) *
  #                            c(1,-1))
  #}
  #coef_g = rep(1.0, len = n_relevant_x)
  #for (ii in 1:n_relevant_x) {
  #  # overlap species
  #  theta[(st:(st + n_relevant_taxa - 1)) %% 4 + 1, 3 * ii - 2] = coef_g[ii] *#al posto di 4 dim
  #    sample(coef)[((ii - 1):(ii + n_relevant_taxa - 2)) %% n_relevant_taxa + 1]
  #  st = st + 1
  #}

  theta[, 1] <- c(-2.0, 1.0, 1.5, 0)
  theta[, 2] <- c(1.0, 2.0, 0, -1.0)
  theta <- t(theta)

  intercept <- matrix(0, n, 4)

  for(i in 1:n){
    if(label[i] == 1){
      intercept[i,] <- mvtnorm::rmvnorm(1, t(x[i,])*(beta1%*%diag(4)), diag(sqrt(0.5), 4))
    }
    if(label[i] == 2){
      intercept[i,] <- mvtnorm::rmvnorm(1, t(x[i,])*(beta2%*%diag(4)), diag(sqrt(0.5), 4))
    }
    if(label[i] == 3){
      intercept[i,] <- mvtnorm::rmvnorm(1, t(x[i,])*(beta3%*%diag(4)), diag(sqrt(0.5), 4))
    }
  }

  Y <- matrix(0, n, 4)

  for(i in 1:n){
    thisrow = as.vector(exp(intercept[i,]+ t(theta)%*% z[i, ]))
    pi = bayess::rdirichlet(n = 1, par = thisrow)
    #pi = thisrow/sum(thisrow)
    #cat("pi", i, ": ", pi, "\n")
    Y[i, ] = rmultinom(1, 1, pi)
  }

  x <- data.frame(x)
  z <- data.frame(z)

  x$X3 <- as.factor(x$X3)
  x$X4 <- as.factor(x$X4)

  #set.seed(1)
  idx <- sample(1:dim(Y)[1], (n*0.75), replace=FALSE)
  Ytrain <- Y[idx,]
  Ytest <- Y[-idx,]
  Xtrain <- x[idx,,drop=FALSE]
  Xtest <- x[-idx,,drop=FALSE]
  Ztrain <- z[idx,,drop=FALSE]
  Ztest <- z[-idx,,drop=FALSE]

  dat$labeltrain <- label[idx]
  dat$labeltest <- label[-idx]
  dat$intercept <- intercept
  dat$intercept_train <- intercept[idx,]
  dat$intercept_test <- intercept[-idx,]
  dat$y <- Y
  dat$Ytrain <- Ytrain
  dat$Ytest <- Ytest
  dat$X <- x
  dat$Xtrain <- Xtrain
  dat$Xtest <- Xtest
  dat$Z <- z
  dat$Ztrain <- Ztrain
  dat$Ztest <- Ztest

  dat$theta <- theta
  dat$nclus <- 3
  return(dat)
}

#' postquant
#'
#' @export
#'
postquant <- function(y, output, data, lab, plot){#, minbinder = F){
  cls <- as.matrix(output$label)
  psm <- mcclust::comp.psm(cls)
  #if(minbinder == T){
  #  mc <- minbinder.ext(psm)
  #} else {
    mc <- minVI(psm)
  #  }
  yhat <- output$pred
  vec <- (c(y)-c(yhat))
  mmse <- mean((y-yhat)^2)
  ari <- adjustedRandIndex(mc$cl, data$label)
  ess <- effectiveSize(output$nclu)

  mypostquant <- list("nclupost" = mean(output$nclu), "MSE" = mmse,
                      "lpml" = output$lpml, "ARI" = ari, "ESS" = ess)
                      #"comp_time" = time[1],
  if(lab==T){
    mypostquant <- list("nclupost" = mean(output$nclu), "MSE" = mmse,
                        "lpml" = output$lpml, "ARI" = ari,
                        #"comp_time" = time[1],
                        "ESS" = ess, "lab" = mc$cl)
  }
  if(plot == T){
    par(mfrow=c(1,2))
    plot(output$nclu, type="l")
    acf(output$nclu)
  }
  return(mypostquant)
}

#' postquant Multinomial Dirichlet
#'
#' @export
#'
postquant_dm <- function(y, yp, output, data, plot, minbinder = F){
  cls <- as.matrix(output$label)
  psm <- comp.psm(cls)
  if(minbinder == T){
    mc <- minbinder.ext(psm)
  } else {
    mc <- minVI(psm)
  }
  ari <- adjustedRandIndex(mc$cl, data$labeltrain)
  ess <- effectiveSize(output$nclu)

  #predictive AUC multiclass
  probs <- apply(output$pipred, c(1, 2), mean)
  colnames(probs) <- colnames(Ytest) <- c(1:4)

  catvec <- c()
  for(i in 1:nrow(Ytest)){
    catvec <- c(catvec, which(yp[i,] == 1))
  }

  auc <- pROC::multiclass.roc(catvec, probs)$auc[1]
  mypostquant <- list("nclupost" = mean(output$nclu),
                      "ARI" = ari, "ESS" = ess, "lab" = mc$cl,
                      "waic" = output$WAIC, "lpml" = output$lpml, "auc" = auc)
  if(plot == T){
    par(mfrow=c(1,2))
    plot(output$nclu, type="l")
    acf(output$nclu)
  }
  return(mypostquant)
}

plot_auc <- function(output_ppm, output_ppmx){
  probs_ppm <- apply(output_ppm$pipred, c(1, 2), mean)
  probs_ppmx <- apply(output_ppmx$pipred, c(1, 2), mean)

  colnames(Ytest) <- c("a_true", "b_true", "c_true", "d_true")
  colnames(probs_ppm) <- c("a_pred_ppm", "b_pred_ppm", "c_pred_ppm", "d_pred_ppm")
  colnames(probs_ppmx) <- c("a_pred_ppmx", "b_pred_ppmx", "c_pred_ppmx", "d_pred_ppmx")
  final_df <- data.frame(cbind(Ytest, probs_ppm, probs_ppmx))

  roc_res <- multi_roc(final_df, force_diag = T)
  plot_roc_df <- plot_roc_data(roc_res)

  ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
    geom_path(aes(color = Group, linetype=Method), size=0.5)
}

####post processing sketch as described in Richardson & Greene (1997)
pp_cs_rg <- function(output, post, nout, dim, refdim = 1){
  C_binder <- max(post$lab)
  median_C <- median(output$nclu)
  nc <- C_binder#median_C#

  eta_work <- matrix(0, nc, dim,)
  eta_pp <- matrix(0, nc, dim)
  noutC <- 0
  for(l in 1:nout){
    if(output$nclu[l]==nc){
      eta_work <- output$eta[c(1:nc),,l]
      eta_pp = eta_pp + eta_work[order(eta_work[,refdim]),]
      noutC <- noutC + 1
    }
  }
  eta_pp <- eta_pp/noutC
  return(eta_pp)
}
