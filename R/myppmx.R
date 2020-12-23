#' PPMx
#' 
#' @param iter MCMC number of iteration
#' @param burn MCMC iteration discarded due to burnin
#' @param thin thinning for MCMC
#' @param nobs length of response vector
#' @param ncon number of continuous covariates
#' @param ncat number of categorical covariates
#' @param catvec \eqn{ncat \times 1} vector indicating the number of categories for each categorical covariate
#' @param alpha value of \eqn{\alpha} for cohesion function (concentration parameter in DP)
#' @param cohesion type of cohesion function to use in the PPMx prior 
#'   1 - Dirichlet process style of cohesion \eqn{c(S) = \alpha \times (\mid S\mid - 1)!}
#'   2 - Uniform cohesion \eqn{c(S) = 1}
#' @param similarity type of similarity function that is employed for the PPMx prior on partitions. Options are
#'   1 - Auxiliary similarity
#'   2 - Double dipper similarity
#' @param consim integer 1 or 2.  1 implies sim for con var is N-N.  2 implies sim is N-NIG
#' @param vector (dimension \eqn{nobs\times 1}) containing response vector
#' @param xcon \eqn{nobs \times ncon} contiguous double vector that contains continuous covariate values
#' @param xcat \eqn{nobs \times ncat} contiguous int vector that contains categorical covariate values
#' @param npred integer indicating number predictions
#' @param xconp  \eqn{npred \times ncon} matrix of continuous covariates to make predictions
#' @param xcatp \eqn{npred \times ncat} matrix of categorical covariates to make predictions
#' @param similparam vector containing similarity functions paramaters
#' @param modelpriors vector containing prior values for model
#' @param mhtune vector containing tuning parameters for MCMC updates
#' @param calibration whether the similarity should be calibrated. Options are
#'   0 - no calibration
#'   1 - standardize similarity value for each covariate
#'   2 - coarsening is applied so that each similarity is raised to the 1/p power
#' @export
# mancano come input tutti gli storage per l output che inizializzo in R e passo come input qui
# poi in R li metto in una lista 
# se non funziona l output di myppmx deve essere una lista

my_ppmx <- function(y,X=NULL,Xpred=NULL,
                          cohesion=1, alpha=1,
                          similarity=1, consim=1,
                          calibration=0,
                          similparam=c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1, 1.0),
                          modelpriors = c(0, 100^2, 1, 1),
                          mhtune=c(0.5, 0.5),
                          iter=1100,burn=100,thin=1){
  
  # X - data.frame whose columns are
  # gcontype - similarity function (1 - Auxilliary, 2 - double dipper)
  # consim=1 => that similarity function of continuous covariate is N-N model (v_j is fixed)
  # consim=2 => that similarity functino of continuous covariate is N-NIG model (v_j is unknown)
  # modelPriors = (mu0, s^2, a1, m)
  # simParms = (m0, s20, v2, k0, nu0, dir, alpha)
  # mh - tuning parameters for MCMC updates of sig2 and sig20
  out <- NULL
  
  #PPM <- ifelse(is.null(X), TRUE, FALSE)
  
  if(!is.data.frame(X) & !is.null(X)) X <- data.frame(X)
  if(!is.data.frame(Xpred) & !is.null(Xpred)){
    Xpred <- data.frame(Xpred);
    colnames(Xpred) <- colnames(X)
  }
  
  nout <- (iter-burn)/thin
  
  nobs <- length(y)
  nxobs <- ifelse(is.null(X), 0, nrow(X))
  npred <- ifelse(is.null(Xpred), 0, nrow(Xpred))
  
  Xall <- rbind(X, Xpred)
  #print(Xall)
  #nmissing <- sum(is.na(Xall))
  #
  #if(nmissing > 0){
  #  Mall <- 1*is.na(Xall) # create the missing matrix with 1 = missing
  #}
  
  # Function that relabels categorical variables to begin with 0
  relab <- function(x) as.numeric(as.factor(as.character(x))) - 1
  
  classes <- sapply(Xall, class)
  catvars <- classes %in% c("factor","character")
  
  # standardize continuous covariates
  if(nxobs > 0){
    if(sum(!catvars) > 0){
      Xconstd <- apply(Xall[,!catvars, drop=FALSE], 2, scale)
      xcon <- Xconstd[1:nobs,,drop=FALSE];
      ncon <- ncol(xcon)
      #if(nmissing > 0) Mcon <- Mall[1:nobs, !catvars, drop=FALSE]
    }else{
      xcon <- cbind(rep(0,nobs));
      ncon <- 0
      #if(nmissing > 0) Mcon <- cbind(rep(0,nobs))
    }
    
    
    if(sum(catvars) > 0){
      # Change the factors or characters into integers with category starting at 0.
      Xcatall <- apply(Xall[, catvars,drop=FALSE], 2, relab)
      #if(nmissing > 0) Mcat <- Mall[1:nobs, catvars, drop=FALSE]
      
      xcat <- Xcatall[1:nobs,,drop=FALSE];
      catvec <- apply(xcat,2,function(x)length(unique(x)))
      ncat <- ncol(xcat)
      
    }else{
      xcat <- cbind(rep(0,nobs));
      #if(nmissing > 0) Mcat <- cbind(rep(0,nobs))
      catvec <- 0
      ncat <- 0
    }
  } else {
    xcon <- cbind(rep(0,1));
    xcat <- cbind(rep(0,1));
    #if(nmissing > 0) Mcon <- cbind(rep(0,1));
    #if(nmissing > 0) Mcat <- cbind(rep(0,1));
  }
  
  
  # Create the matrices of continuous and categorical variables
  # that will be used in th posterior predictive distributions
  if(npred > 0){
    if(sum(!catvars) > 0){
      Xconstd <- apply(Xall[,!catvars, drop=FALSE], 2, scale)
      xconp <- Xconstd[(nrow(Xall)-npred+1):nrow(Xall),,drop=FALSE];
      #if(nmissing > 0) Mconp <- Mall[(nrow(Xall)-npred+1):nrow(Xall), !catvars, drop=FALSE]
      ncon <- ncol(xconp)
    } else {
      xconp <- cbind(rep(0,npred));
      #if(nmissing > 0) Mconp <- cbind(rep(0,npred));
      ncon <- 0
    }
    
    if(sum(catvars) > 0){
      Xcatall <- apply(Xall[, catvars,drop=FALSE], 2, relab)
      xcatp <- Xcatall[(nrow(Xall)-npred+1):nrow(Xall),,drop=FALSE];
      #if(nmissing > 0) Mcatp <- Mall[(nrow(Xall)-npred+1):nrow(Xall), catvars, drop=FALSE]
      ncat <- ncol(xcatp)
      catvec <- apply(Xcatall,2,function(x)length(unique(x)))
    } else {
      xcatp <- cbind(rep(0,npred));
      ncat <- 0
      catvec <- 0
      #if(nmissing > 0) Mcatp <- cbind(rep(0,npred));
    }
    
  } else {
    xconp <- cbind(rep(0,1));
    xcatp <- cbind(rep(0,1));
    #if(nmissing > 0) Mconp <- cbind(rep(0,1));
    #if(nmissing > 0) Mcatp <- cbind(rep(0,1));
  }
  
  dissimtn <- 0
  dissimtt <- 0
  
  # Create empty vectors that will hold MCMC iterates
  mu <- sig2 <- Si <- like <- ispred <- zi <- isordpred <- matrix(1,nrow=nout,ncol=nobs)
  mu0 <- sig20 <- nclus <- rep(1,nout)
  ppred <- predclass <-  rbpred <-  matrix(1, nrow=nout, ncol=npred)
  predclass_prob <- matrix(1, nrow=nout, ncol=npred*nobs)
  WAIC <- lpml <- rep(1,1)
  
  out <- myppmx(as.integer(iter), as.integer(burn), as.integer(thin), 
                  as.integer(nobs), as.integer(ncon), as.integer(ncat), 
                  as.vector(catvec), as.double(alpha), as.integer(cohesion), 
                  as.integer(similarity), as.integer(consim), as.vector(y), as.vector(xcon), 
                  as.vector(xcat), as.integer(npred), as.matrix(xconp), as.matrix(xcatp), 
                  as.vector(similparam), as.vector(modelpriors), as.vector(mhtune), 
                  as.integer(calibration))
  
  res <- list()
  res$mu <- matrix(out$mu, nrow=nout, byrow=TRUE)
  res$sig2 <- matrix(out$sig2, nrow=nout, byrow=TRUE)
  res$Si <- matrix(out$Si, nrow=nout, byrow=TRUE)
  res$like <- matrix(out$like, nrow=nout, byrow=TRUE)
  res$fitted <- matrix(out$ispred, nrow=nout, byrow=TRUE)
  res$ppred <- matrix(out$ppred, nrow=nout, byrow=TRUE)
  res$predclass <- matrix(out$predclass, nrow=nout, byrow=TRUE)
  res$predclass_prob <- matrix(out$predclass_prob, nrow=nout, byrow=TRUE)
  res$rbpred <- matrix(out$rbpred, nrow=nout, byrow=TRUE)
  res$mu0 <- out$mu0
  res$sig20 <- out$sig20
  res$nclus <- out$nclus
  res$WAIC <- out$WAIC
  res$lpml <- out$lpml
  
  return(res)
}

