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
  out$nj <- nj[1:nclus]
  return(out)
}

genera_dati <- function(n = 100, P = 2, Q = 2, dim =2){
  discretize <- function(vec){
    vec[which(vec<0)] <- 0
    vec[which(vec>0)] <- 1
    #vec <- as.factor(vec)
    return(vec)
  }
  ncov <- Q + P
  ##prima ci stanno le binarie
  XX <- matrix(0, n, (ncov))

  Y <- matrix(0, n, dim)
  beta1 <- matrix(1, ncov, dim)
  beta2 <- matrix(-1.3, ncov, dim)
  beta3 <- matrix(1.3, ncov, dim)

  pro <- c(0.2,0.5,0.3)
  clusterlabel <- c()
  for(i in 1:n){
    u <- runif(1)
    if(u<pro[1]){
      #cat(i,"ciao1","\n")
      #X e Z li devo costruire qui
      if(Q!=0){XX[i, c(1:Q)] <- discretize(rmvnorm(1, rnorm(Q, -2.1, .25)))}
      if(P!=0){XX[i, c((Q+1):ncov)] <- rmvnorm(1, rnorm(P, -2.1, .25))}
      #Y[i,] <- mvtnorm::rmvnorm(1, mean=XX[i,]%*%beta1, sigma = diag(.1, nrow=dim))
      #clusterlabel[i] <- 1
    }
    else{
      if(u<(pro[1]+pro[2])){
        #cat(i,"ciao2","\n")
        if(Q!=0){XX[i, c(1:Q)] <- discretize(rmvnorm(1, rnorm(Q, 0.0, .25)))}
        if(P!=0){XX[i, c((Q+1):ncov)] <- rmvnorm(1, rnorm(P, 0.0, .25))}
        #Y[i,] <- mvtnorm::rmvnorm(1, mean=XX[i,]%*%beta1, sigma = diag(.1, nrow=dim))
        #clusterlabel[i] <- 2
      }
      else{
        #cat(i,"ciao3","\n")
        if(Q!=0){XX[i, c(1:Q)] <- discretize(rmvnorm(1, rnorm(Q, 2.1, .25)))}
        if(P!=0){XX[i, c((Q+1):ncov)] <- rmvnorm(1, rnorm(P, 2.1, .25))}
        #Y[i,] <- mvtnorm::rmvnorm(1, mean=XX[i,]%*%beta1, sigma = diag(.1, nrow=dim))
        #clusterlabel[i] <- 3
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
gcd <- function(n_obs = 100, concov = 2, K = 2, similarity = 1, simparm = 1,
                alpha = 1, m0 = 0, s20 = 1, v = 2, k0 = 10, v0 = 1, plot = F){
  myinfopart <- NULL
  concov = 2
  n_obs = 100
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

#' postquant
#'
#' @export
#'

postquant <- function(y, output, data, lab, plot, minbinder = F){
  cls <- as.matrix(output$label)
  psm <- comp.psm(cls)
  if(minbinder == T){
    mc <- minbinder.ext(psm)
  } else {
    mc <- minVI(psm)
    }
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
