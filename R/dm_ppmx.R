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
#' @param CC number of auxiliary parameters
#' @param reuse option for the reuse algorithm by Favaro&Teh. integer 0 or 1.
#'   0 - reuse algorithm is NOT adopted
#'   1 - reuse algorithm is adopted
#' @param similarity type of similarity function that is employed for the PPMx prior on partitions. Options are
#'   1 - Auxiliary similarity
#'   2 - Double dipper similarity
#' @param consim integer 1 or 2.  1 implies sim for con var is N-N.  2 implies sim is N-NIG
#' @param vector (dimension \eqn{nobs\times 1}) containing response vector
#' @param xcon \eqn{nobs \times ncon} contiguous double vector that contains continuous covariate values
#' @param xcat \eqn{nobs \times ncat} contiguous int vector that contains categorical covariate values
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

my_dm_ppmx <- function(y, X=NULL, alpha=1, CC = 3, reuse = 1, PPMx = 1, similarity = 1, consim=1, calibration=0,
                        similparam=c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1, 1.0),
                        modelpriors,
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


  if(!is.data.frame(X) & !is.null(X)) X <- data.frame(X)

  nout <- (iter-burn)/thin

  nobs <- dim(y)[1]
  nxobs <- ifelse(is.null(X), 0, nrow(X))

  Xall <- rbind(X)

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
      catvec <- 0
      ncat <- 0
    }
  } else {
    xcon <- cbind(rep(0,1));
    xcat <- cbind(rep(0,1));
  }

  alpha <- alpha#similparam[7]
  hP0_m0 <- as.vector(modelpriors$hP0_m0)
  hP0_L0 <- as.vector(modelpriors$hP0_L0)
  hP0_nu0 <- as.double(modelpriors$hP0_nu0)
  hP0_V0 <- as.vector(modelpriors$hP0_V0)

  out <- dm_ppmx(as.integer(iter), as.integer(burn), as.integer(thin),
                  as.integer(nobs), as.integer(PPMx), as.integer(ncon), as.integer(ncat),
                  as.vector(catvec), as.double(alpha), as.integer(CC), as.integer(reuse),
                  as.integer(consim), as.integer(similarity),
                  as.integer(calibration), as.matrix(y),
                  as.vector(xcon), as.vector(xcat), as.vector(similparam),
                  as.vector(hP0_m0), as.vector(hP0_L0), as.double(hP0_nu0),
                  as.vector(hP0_V0), as.vector(mhtune))

  res <- NULL
  nclu <- out$nclus
  res$nclu <- nclu

  nclu_cs <- cumsum(nclu)
  #mu_out <- out$mu#matrix(out$mu, nrow=nout*nobs, byrow=TRUE)
  #mu_ar <- array(0, dim = c(max(nclu), ncol(y), nout))
  #for(l in 1:nout){
  #  for(i in 1:nclu[l]){
  #    mu_ar[i, ,1] <- mu_out[i, ]
  #    if(l > 1){
  #      mu_ar[i, ,l] <- mu_out[i+nclu_cs[l-1], ]
  #    }
  #  }
  #}

  eta_out <- out$eta#matrix(out$mu, nrow=nout*nobs, byrow=TRUE)
  eta_ar <- array(0, dim = c(max(nclu), ncol(y), nout))
  for(l in 1:nout){
    for(i in 1:nclu[l]){
      eta_ar[i, ,1] <- eta_out[i, ]
      if(l > 1){
        eta_ar[i, ,l] <- eta_out[i+nclu_cs[l-1], ]
      }
    }
  }

  pi_out <- out$pi

  res$label <- matrix(out$cl_lab, nrow = nout, byrow=TRUE)
  #res$mu <- mu_ar
  res$eta <- eta_ar
  ypred <- array(0, dim = c(nrow(y), ncol(y), nout))

  for(l in 1:nout){
    for(i in 1:nrow(y)){
      ypred[i,,l] <- rmultinom(1, 1, pi_out[i, ,l])
    }
    #caret::confusionMatrix(as.factor(ypred),as.factor(y))
  }

  #ypred <- apply(ypred, c(1, 2), mean)
  #res$pred <- t(apply(out$ispred, c(1, 2), mean))
  #res$fitted <- matrix(out$ispred, nrow=nout, byrow=TRUE)
  #res$nclus <- out$nclus
  #res$WAIC <- out$WAIC
  res$lpml <- out$lpml

  return(res)
}

#likelihood e lpml in cpp

#funzione x gen dati
#script
