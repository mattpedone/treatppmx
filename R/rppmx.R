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
  
  cat("xcon", as.vector(t(xcon)), "\n")
  cat("xcat", as.vector(t(xcat)), "\n")
  cat("Cvec", as.vector(t(Cvec)), "\n")
  cat("m0", as.double(m0), "\n")
  cat("k0", as.double(k0), "\n")
  cat("v0", as.double(v0), "\n")
  cat("s20", as.double(s20), "\n")
  cat("v", as.double(v), "\n")
  cat("dw", as.vector(dirweights), "\n")
  
  Cout <- rppmx(as.integer(nobs), as.integer(similarity), as.integer(simparm), as.double(alpha), 
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