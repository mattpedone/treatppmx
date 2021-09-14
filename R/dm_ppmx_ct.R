#' PPMx
#'
#' @param y response variable
#' @param X covariates
#' @param Xpred covariate for prediction
#' @param z prognostic covariates
#' @param zpred prognostic covariates for prediction
#' @param asstreat treatment for patients in training set. vector of charachters
#' @param iter MCMC number of iteration
#' @param burn MCMC iteration discarded due to burnin
#' @param thin thinning for MCMC
#' @param nobs length of response vector
#' @param ncon number of continuous covariates
#' @param ncat number of categorical covariates
#' @param catvec \eqn{ncat \times 1} vector indicating the number of categories for each categorical covariate
#' @param alpha value of \eqn{\alpha} for cohesion function (concentration parameter in DP)
#' @param sigma value of \eqn{\sigma} for cohesion function (reinforcement parameter in NGG)
#' @param CC number of auxiliary parameters
#' @param reuse option for the reuse algorithm by Favaro&Teh. integer 0 or 1.
#'   0 - reuse algorithm is NOT adopted
#'   1 - reuse algorithm is adopted
#' @param cohesion type of cohesion function that is employed for the PPMx prior on partitions. Options are
#'   1 - DirichletProcess-like cohesion (DP) cohesion
#'   2 - Normalized Generalized Gamma Process (NGG) cohesion
#' @param similarity type of similarity function that is employed for the PPMx prior on partitions. Options are
#'   1 - Auxiliary similarity
#'   2 - Double dipper similarity
#'   3 - Gower dissimilarity
#' @param gowtot if similarity parameter is 3 then Gower dissimilarity is employed. default total.
#' if gowtot ==0,  then Average Gower dissimilarity is taken
#' @param alphagow $\alpha$ parameter to compute gower dissimilarity
#' @param consim integer 1 or 2.  1 implies sim for con var is N-N.  2 implies sim is N-NIG
#' @param vector (dimension \eqn{nobs\times 1}) containing response vector
#' @param xcon \eqn{nobs \times ncon} contiguous double vector that contains continuous covariate values
#' @param xcat \eqn{nobs \times ncat} contiguous int vector that contains categorical covariate values
#' @param similparam vector containing similarity functions paramaters
# double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
#' @param modelpriors vector containing prior values for model
# @param mhtune vector containing tuning parameters for MCMC updates
#' @param calibration whether the similarity should be calibrated. Options are
#'   0 - no calibration
#'   1 - standardize similarity value for each covariate
#'   2 - coarsening is applied so that each similarity is raised to the 1/p power
#' @param coardegree option to temper the coarsening
#'   1 - $g(x^\star)^{1/p}$
#'   2 - $g(x^\star)^{1/p^{1/2}}$
#'   3 - $g(x^\star)^{1/p^{1/3}}$
#' @param update_hierarchy should hypeparameter for BNP intercept be updated? if 1 yes (default)
#' @param initpc: partial correlation initialization for prognostic covariates coefficient (default)
#' @export
# mancano come input tutti gli storage per l output che inizializzo in R e passo come input qui
# poi in R li metto in una lista
# se non funziona l output di myppmx deve essere una lista

my_dm_ppmx_ct <- function(y, X=NULL, Xpred = NULL, z=NULL, zpred=NULL, asstreat = NULL, alpha=1.0,
                       sigma = 0.2, CC = 3, reuse = 1, PPMx = 1, cohesion = 2, similarity = 1, consim=1,
                       gowtot = 1, alphagow = 1, calibration=0, coardegree = 1,
                       similparam, modelpriors, update_hierarchy = 1, hsp = 1, iter=1100,
                       burn=100,thin=1, mhtunepar = c(.05, .05), nclu_init = 5){

  # X - data.frame whose columns are
  # gcontype - similarity function (1 - Auxilliary, 2 - double dipper)
  # consim=1 => that similarity function of continuous covariate is N-N model (v_j is fixed)
  # consim=2 => that similarity functino of continuous covariate is N-NIG model (v_j is unknown)
  # modelPriors = (mu0, s^2, a1, m)
  # simParms = (m0, s20, v2, k0, nu0, dir, alpha)
  # mh - tuning parameters for MCMC updates of sig2 and sig20


  out <- NULL

  if(!is.data.frame(X) & !is.null(X)) X <- data.frame(X)
  if(!is.data.frame(Xpred) & !is.null(Xpred)){
    Xpred <- data.frame(Xpred);
    colnames(Xpred) <- colnames(X)
  }
  cnames <- colnames(X)

  nout <- (iter-burn)/thin
  nobs <- dim(y)[1]

  if(is.null(X) & is.null(Xpred)){
    xcon <- cbind(rep(0,1));
    xconp <- cbind(rep(0,1));
    xcat <- cbind(rep(0,1));
    xcatp <- cbind(rep(0,1));
    ncon <- 0
    ncat <- 0
    npred <- 0
  }

  if(!(is.null(X) & is.null(Xpred))){
    nxobs <- ifelse(is.null(X), 0, nrow(X))
    npred <- ifelse(is.null(Xpred), 0, nrow(Xpred))
    Xall <- rbind(X, Xpred)
    # Function that relabels categorical variables to begin with 0
    relab <- function(x) as.numeric(as.factor(as.character(x))) - 1
    classes <- sapply(Xall, class)
    catvars <- classes %in% c("factor","character")

    # standardize continuous covariates
    if(nxobs > 0){
      if(sum(!catvars) > 0){
        Xconstd <- Xall[,!catvars, drop=FALSE]#apply(Xall[,!catvars, drop=FALSE], 2, scale)#
        xcon <- Xconstd[1:nobs,,drop=FALSE];
        #print(xcon)
        ncon <- ncol(xcon)
      }else{
        xcon <- cbind(rep(0,nobs));
        ncon <- 0
      }

      if(sum(catvars) > 0){
        # Change the factors or characters into integers with category starting at 0.
        Xcatall <- apply(Xall[, catvars,drop=FALSE], 2, relab)
        xcat <- Xcatall[1:nobs,,drop=FALSE];
        catvec <- apply(xcat,2,function(x)length(unique(x)))
        ncat <- ncol(xcat)
      }else{
        xcat <- cbind(rep(0,nobs));
        catvec <- 0
        ncat <- 0
      }
    }
    # Now consider the case when number of covariates for prediction are greater than zero
    if(npred > 0){
      if(sum(!catvars) > 0){
        Xconstd <- Xall[,!catvars, drop=FALSE]#apply(Xall[,!catvars, drop=FALSE], 2, scale)#
        xconp <- Xconstd[(nrow(Xall)-npred+1):nrow(Xall),,drop=FALSE];
        #print(xconp)
        ncon <- ncol(xconp)
      } else {
        xconp <- cbind(rep(0,npred));
        ncon <- 0
      }

      if(sum(catvars) > 0){
        Xcatall <- apply(Xall[, catvars,drop=FALSE], 2, relab)
        xcatp <- Xcatall[(nrow(Xall)-npred+1):nrow(Xall),,drop=FALSE];
        ncat <- ncol(xcatp)
        catvec <- apply(Xcatall,2,function(x)length(unique(x)))
      } else {
        xcatp <- cbind(rep(0,npred));
        ncat <- 0
        catvec <- 0
      }
    }
  }

  if(similarity == 3){
    dissim <- as.matrix(cluster::daisy(Xall, metric="gower"))
    dissimtn <- dissim[1:nobs, 1:nobs]
    dissimtt <- dissim[-c(1:nobs), 1:nobs]
  }else{
    dissimtn <- 0
    dissimtt <- 0
  }

  cormat = matrix(0, ncol(y), ncol(z))
  pmat = matrix(0, ncol(y), ncol(z))
  yy = y/rowSums(y) # compositionalize
  for(rr in 1:ncol(y)){
    for(cc in 1:ncol(z)){
      pmat[rr, cc] = suppressWarnings(stats::cor.test(z[, cc], yy[, rr], method = "spearman",
                                     exact = F)$p.value)
      cormat[rr, cc] = suppressWarnings(stats::cor(z[, cc], yy[, rr], method = "spearman"))
    }
  }
  # defaults to 0.2 false discovery rate
  pm = matrix((stats::p.adjust(c(pmat), method = "fdr") <= 0.2) + 0, ncol(y),
              ncol(z))
  betmat = cormat * pm
  beta <- as.vector(t(betmat))
  beta[is.na(beta)] <- 0.0
  beta <- beta + 0.0001

  A <- length(unique(asstreat))#nT int
  treatments <- asstreat#treatments vec
  n_a <- table(treatments)#num_treat vec 2
  max_n_treat <- max(n_a)#num_treat.max

  curr_cluster <- matrix(0, A, max_n_treat)
  card_cluster <- matrix(0, A, max_n_treat)
  nclu_curr <- rep(nclu_init, A)
  for(a in 1:A){
    work <- kmeans(X[treatments==(a),], nclu_init, iter.max = 10, nstart = 25)
    for(i in 1:length(work$cluster)){
      curr_cluster[a,i] <- work$cluster[i]
    }
    for(i in 1:nclu_init){
      card_cluster[a,i] <- table(work$cluster)[i]
    }
  }

  treatments <- treatments - 1

  #calcola 1 colonna V
  #calcola resto matrice
  #passa a funzione cpp
  #   -V
  #   -sigma
  #   -cohesion

  a <- alpha
  sigma <- sigma
  integrand <- function(u, n, sigma, a){u^(n-1)*exp((a/sigma)-((a*((1+u)^sigma))/(sigma)))*((1+u)^(-n+sigma))}
  Vwm <- matrix(NA, nrow = max_n_treat+2, ncol = max_n_treat +2)
  Vwm[1, 1] <- 1

  for(n in 2:(max_n_treat+2)){
    numsol <- as.double(integrate(integrand, lower = 0, upper = Inf, n = n, sigma = sigma, a = a)[1])
    Vwm[n, 1] <- (a/gamma(n))*numsol
    }

  for(n in 1:max_n_treat+1){
    for(k in 1:n){
      Vwm[n+1, k+1] <- Vwm[n, k] -((n-(sigma*k))*Vwm[n+1, k])
    }
  }

  alpha <- alpha#similparam[7]
  hP0_m0 <- as.vector(modelpriors$hP0_m0)
  hP0_L0 <- as.vector(modelpriors$hP0_L0)
  hP0_nu0 <- as.double(modelpriors$hP0_nu0)
  hP0_V0 <- as.vector(modelpriors$hP0_V0)
  treatments <- c(as.vector(treatments))

  out <- dm_ppmx_ct(as.integer(iter), as.integer(burn), as.integer(thin),
                 as.integer(nobs), as.vector(treatments), as.integer(PPMx),
                 as.integer(ncon), as.integer(ncat),
                 as.vector(catvec), as.double(alpha), as.double(sigma), as.matrix(Vwm), as.integer(cohesion),
                 as.integer(CC), as.integer(reuse),
                 as.integer(consim), as.integer(similarity), as.integer(gowtot),
                 as.integer(alphagow), as.vector((dissimtn)), as.vector((dissimtt)),
                 as.integer(calibration), as.integer(coardegree), as.matrix(y), as.matrix(z), as.matrix(zpred),
                 as.vector(t(xcon)), as.vector(t(xcat)), as.vector(t(xconp)),
                 as.vector(t(xcatp)), as.integer(npred), as.vector(similparam),
                 as.vector(hP0_m0), as.vector(hP0_L0), as.double(hP0_nu0),
                 as.vector(hP0_V0), as.integer(update_hierarchy), as.vector(t(beta)),
                 as.integer(hsp), as.vector(mhtunepar), as.integer(A), as.vector(n_a),
                 as.matrix(curr_cluster), as.matrix(card_cluster), as.vector(nclu_curr))

  ###PREPARE OUTPUT
  res <- list()

  #label (clustering)
  label <- list()
  for(t in 1:A){
    label[[t]] <- out$cl_lab[t,,]
  }

  #names(label) <- treatnames
  res$label <- label
  res$asstreat <- asstreat
  num_treat <- out$num_treat
  #names(num_treat) <- treatnames
  res$num_treat <- num_treat

  #number of cluster
  nclu <- out$nclus
  res$nclu <- nclu#(nT x nout matrix)

  nclu_cs <- apply(nclu, 1, cumsum)#(nt x nout matrix)


  #intercept (and its acceptance rate)
  eta_out <- out$eta#matrix(out$mu, nrow=nout*nobs, byrow=TRUE)
  #print(eta_out)
  #eta_ar <- array(0, dim = c(max(nclu), ncol(y), nout, nT))
  #for(t in 1:nT){
  #  for(l in 1:nout){
  #    for(i in 1:nclu[l, t]){
  #      eta_ar[i, ,1, t] <- eta_out[i, ,t]
  #      if(l > 1){
  #        eta_ar[i, ,l, t] <- eta_out[i+nclu_cs[l-1], , t]
  #      }
  #    }
  #  }
  #}
  res$eta <- eta_out
  #res$acc_rate_eta <- eta$acc#sum(out$eta_acc)#/sum(out$nclu)

  #prognostic covariates
  #print((Ytest))
  beta <- array(0, dim = c(ncol(z), length(y), nout))

  beta0 <- out$beta
  for(iter in 1:nout){
    for(k in 1:ncol(Y)){
      for(q in 1:ncol(z)){
        h <- q + (k-1) * ncol(z)
        beta[q, k, iter] <- beta0[h, iter]
      }
    }
  }

  res$beta <- beta
  res$acc_rate_beta <- out$beta_acc/iter#out$beta_acc#/sum(out$nclu)
  #pi (dirichlet parameter)
  pi_out <- out$pi
  res$pi_out <- pi_out
  res$pipred <- out$pipred

  #in sample prediction & model fit
  #yispred <- array(0, dim = c(nobs, ncol(y), nout))
  res$isypred <- out$yispred
  res$WAIC <- out$WAIC
  res$lpml <- out$lpml

  #posterior predictive
  ypred <- array(0, dim = c(npred, ncol(y), A, nout))
  for(i in 1:nout){
    for(a in 1:A){
      ypred[,,a,i] <- out$ypred[i,1][[1]][,,a]
    }
  }
  res$ypred <- ypred#out$ypred

  #cluster classification prediction (unsupervised clustering)
  clupred <- matrix(0, nout, npred)

  ll = 0
  for(l in 1:nout){
    for(pp in 1:npred){
      clupred[l, pp] <- out$clupred[ll*npred + pp]
      ll = ll+1
    }
  }

  res$clupred <- clupred

  return(res)
}

#likelihood e lpml in cpp

#funzione x gen dati
#script
