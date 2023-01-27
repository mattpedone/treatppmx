#' ppmxct_fix
#'
#' Function to predict personalized treatment for \eqn{n_{test}} new untreated patients, given their biomarkers.
#' It leverages response to treatment, \eqn{P} prognostic, \eqn{Q} predictive biomarkers
#' of \eqn{n_{train}} historical patients.
#' It accounts for \eqn{K} ordinal response level and \eqn{T} competing treatments.
# For any methodological detail refer to Pedone M., Argiento R., and Stingo F.
#'
#' @param y \eqn{n_{train} \times K} matrix of ordinal-valued response variable
#' @param X \eqn{n_{train} \times Q} dataframe of predictive biomarkers
#' @param Xpred \eqn{n_{test} \times Q} dataframe of predictive covariates of new untreated patient
#' @param Z \eqn{n_{train} \times P} dataframe of prognostic covariates
#' @param Zpred \eqn{n_{test} \times P} dataframe of prognostic covariates of new untreated patient
#' @param asstreat \eqn{n_test} vector of integers encoding treatment received by historical patients
#' @param PPMx logical. option for the use of product partition model with covariates. (default is yes)
#' @param cohesion type of cohesion function that is employed for the PPMx prior on partitions. Options are
#'   1 - DirichletProcess-like cohesion (DP) cohesion
#'   2 - Normalized Generalized Gamma Process (NGG) cohesion
#' @param kappa vector of possible values for \eqn{\kappa} for cohesion function (concentration parameter in DP and NGG)
#' @param sigma vector of possible value for \eqn{\sigma} parameter in the cohesion function (reinforcement parameter in NGG)
#' @param similarity type of similarity function that is employed for the PPMx prior on partitions. Options are
#'   1 - Auxiliary similarity
#'   2 - Double dipper similarity
#   3 - Gower dissimilarity
#' @param consim integer 1 or 2.  1 implies sim for con var is NN.  2 implies sim is NNIG
#' @param similparam vector containing similarity functions paramaters
# double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
#' @param calibration If the similarity function is Auxiliary or Double Dipper, the similarity
#' can be calibrated. Options are
#'   0 - no calibration
#'   1 - standardize similarity value for each covariate
#'   2 - coarsening is applied so that each similarity is raised to the 1/p power
#' @param coardegree If the similarity is coarsened, it is possible to temper the coarsening
#'   1 - \eqn{g(x^*)^{1/p}}
#'   2 - \eqn{g(x^*)^{1/p^{1/2}}}
# @param gowtot if similarity parameter is 3 then Gower dissimilarity is employed. default total.
# if gowtot == 0,  then Average Gower dissimilarity is taken
# @param alphagow \eqn{\alpha} parameter to compute gower dissimilarity
#' @param modelpriors vector containing prior values for model
#' @param update_hierarchy should hyperparameter for BNP intercept be updated? if 1 yes (default)
#' @param hsp parameter for employ horseshoe prior for coefficients for prognostic markers
#' @param iter MCMC number of iteration
#' @param burn MCMC iteration discarded due to burnin
#' @param thin thinning for MCMC
#' @param mhtunepar vector containing tuning parameters for MCMC updates
#' @param CC number of auxiliary parameters for Algorithm 8 by Neal (2000)
#' @param reuse option for the reuse algorithm by Favaro and Teh (2013). integer 0 or 1.
#'   0 - reuse algorithm is not adopted
#'   1 - reuse algorithm is adopted
#' @param fix_clu fixed clustering
#'
#' @references
#'
#' Page, G. L. and Quintana, F. A. (2016). Calibrating covariate informed
#' product partition models. \emph{Statistics and Computing}, \strong{28}(5), 1009–1031.
#'
#' Neal, R. M. (2000). Markov chain sampling methods for Dirichlet process mixture
#' models. \emph{Journal of Computational and Graphical Statistics}, \strong{9}(2): 249–265.
#'
#' Favaro, S., Teh, Y. W., et al. (2013). MCMC for normalized random measure mixture
#' models. \emph{Statistical Science}, \strong{28}(3): 335–359.

#' @return a list of 16 elements.
#'   \itemize{
#'   \item \code{label}: List of \eqn{T} matrices. Each element is a \eqn{n_{test}^{a}\times}\code{nout} matrix.
#'   \code{nout} is the number of MCMC iterations after burnin and thinning
#'   The row of each matrix contains the vector of cluster labels of the historical
#'   patient assigned to given treatment
#'   \item \code{asstreat}: Vector of \eqn{n_{test}} treatment assignment
#'   \item \code{num_treat}: Vector of dimension \eqn{T}. Contains the number
#'   of patients assigned to each treatment \eqn{n^a}, for \eqn{a=1, \dots, T}.
#'   \item \code{nclu}: \eqn{T\times}\code{nout} matrix of total number of cluster at each MCMC iteration for each treatment
#'   \item \code{eta}: \code{nout}\eqn{\times T} List of matrices. Each matrix has dimension \eqn{K\times n^a}.
#'   Each element of the list is the linear predictor matrix for treatment \eqn{a} at each MCMC iteration.
#'   \item \code{beta}: Array of dimensions \eqn{P\times K \times}\code{nout}. It contains the matrix of prognostic factors at each MCMC iteration
#'   \item \code{acc_beta}: \eqn{P\times K} matrix containing the counts each prognostic coefficient has been accepted in the MH step. Note that these counts should be divided by \code{iter}
#'   Each element of the list is a list of 2 matrices (one for each treatment).
#'   \item \code{pi_out}: \eqn{n\times K \times}\code{nout} array. It contains the Multinomial parameters for each patients at each MCMC iteration
#'   \item \code{isypred}: \eqn{n\times K \times}\code{nout} array. It contains the outcome in-sample-prediction for each patients at each MCMC iteration
#'   \item \code{WAIC}: scalar. It is the average WAIC
#'   \item \code{lpml}: scalar. It is the average lpml
#'   \item \code{sigmangg}: \eqn{T\times}\code{nout} matrix. It contains the \eqn{\sigma} parameter of NGGP for each treatment at each MCMC iteration
#'   \item \code{kappangg}: \eqn{T\times}\code{nout} matrix. It contains the \eqn{\kappa} parameter of NGGP for each treatment at each MCMC iteration
#'   \item \code{ypred}: \eqn{n_{pred} \times K \times T \times}\code{nout} array. It contains the predicted outcome for each untreated patients, under all the competing for each treatment at each MCMC iteration
#'   \item \code{ypred}: \eqn{n_{pred} \times K \times T \times}\code{nout} array. It contains the Multinomial parameters for each untreated patients, under all the competing for each treatment at each MCMC iteration
#' }
#' @export

ppmxct_fixed <- function(y, X=NULL, Xpred = NULL, Z=NULL, Zpred=NULL, asstreat = NULL,
                   PPMx = 1, cohesion = 2, kappa = c(1.0, 30.0, 10, 1),
                   sigma = c(0.005, 1.0, 10), similarity = 1, consim=1,
                   similparam, calibration=0, coardegree = 1, modelpriors,
                   update_hierarchy = 1, hsp = 1, iter=1100, burn=100, thin=1,
                   mhtunepar = c(.05, .05), CC = 3, reuse = 1, fix_clu){

  # X - data.frame whose columns are
  # gcontype - similarity function (1 - Auxilliary, 2 - double dipper)
  # consim=1 => that similarity function of continuous covariate is N-N model (v_j is fixed)
  # consim=2 => that similarity functino of continuous covariate is N-NIG model (v_j is unknown)
  # modelPriors = (mu0, s^2, a1, m)
  # simParms = (m0, s20, v2, k0, nu0, dir, kappa)
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

    if(nxobs > 0){
      if(sum(!catvars) > 0){
        Xconstd <- Xall[,!catvars, drop=FALSE]
        xcon <- Xconstd[1:nobs,,drop=FALSE]
        #print(xcon)
        ncon <- ncol(xcon)
      }else{
        xcon <- cbind(rep(0,nobs))
        ncon <- 0
      }

      if(sum(catvars) > 0){
        # Change the factors or characters into integers with category starting at 0.
        Xcatall <- apply(Xall[, catvars,drop=FALSE], 2, relab)
        xcat <- Xcatall[1:nobs,,drop=FALSE];
        catvec <- apply(xcat,2,function(x)length(unique(x)))
        ncat <- ncol(xcat)
      }else{
        xcat <- cbind(rep(0,nobs))
        catvec <- 0
        ncat <- 0
      }
    }
    # Now consider the case when number of covariates for prediction are greater than zero
    if(npred > 0){
      if(sum(!catvars) > 0){
        Xconstd <- Xall[,!catvars, drop=FALSE]
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

  pmat = matrix(0, ncol(y), ncol(Z))
  noprog <- 0
  if((is.null(Z)) & (is.null(Zpred))){
    beta <- pmat
    Z <- matrix(0, nrow = nrow(X), ncol = 1)
    Zpred <- matrix(0, nrow = nrow(Xpred), ncol = 1)
    noprog <- 1
  } else {
    cormat = matrix(0, ncol(y), ncol(Z))
    yy = y/rowSums(y) # compositionalize
    for(rr in 1:ncol(y)){
      for(cc in 1:ncol(Z)){
        pmat[rr, cc] = suppressWarnings(stats::cor.test(Z[, cc], yy[, rr], method = "spearman",
                                                        exact = F)$p.value)
        cormat[rr, cc] = suppressWarnings(stats::cor(Z[, cc], yy[, rr], method = "spearman"))
      }
    }
    # defaults to 0.2 false discovery rate
    pm = matrix((stats::p.adjust(c(pmat), method = "fdr") <= 0.2) + 0, ncol(y),
                ncol(Z))
    betmat = cormat * pm
    beta <- as.vector(t(betmat))
    beta[is.na(beta)] <- 0.0
    beta <- beta + 0.0001
  }

  A <- length(unique(asstreat))#nT int
  treatments <- asstreat#treatments vec
  n_a <- table(treatments)#num_treat vec 2
  max_n_treat <- max(n_a)#num_treat.max
  min_n_treat <- min(n_a)#num_treat.max

  #curr_cluster <- matrix(0, A, max_n_treat)
  #card_cluster <- matrix(0, A, max_n_treat)
  #nclu_curr <- rep(nclu_init, A)
  #for(a in 1:A){
  #  work <- kmeans(X[treatments==(a),], nclu_init, iter.max = 10, nstart = 25)
  #  for(i in 1:length(work$cluster)){
  #    curr_cluster[a,i] <- work$cluster[i]
  #  }
  #  for(i in 1:nclu_init){
  #    card_cluster[a,i] <- table(work$cluster)[i]
  #  }
  #}

  #load(clu_fix)

  curr_cluster <- t(fix_clu)
  nclu_curr <- apply(curr_cluster, 1, max)
  card_cluster <- matrix(0, A, max(nclu_curr))
  for(a in 1:A){
    tabcard <- table(curr_cluster[a,])
    for(i in 1:nclu_curr[a]){
      card_cluster[a,i] <- tabcard[i]
    }
  }

  treatments <- treatments - 1

  kappadp <- kappa[4]
  nsigma <- sigma[3]
  nkappa <- kappa[3]
  sigmagrid <- seq(sigma[1], sigma[2], length.out = nsigma)
  kappagrid <- seq(kappa[1], kappa[2], length.out = nkappa)
  grid <- t(expand.grid(kappagrid, sigmagrid))[c(2,1),]
  ngrid <- ncol(grid)

  if(cohesion == 2){
    Vwm <- array(NA, dim = c(max_n_treat+1, max_n_treat+1, ngrid))
    Vwm[1, 1,] <- 1
    for(l in 1:ngrid){
      for(n in (min_n_treat-1):(max_n_treat+1)){
        Vwm[n, (1:n), l] <- vweights::computev(n, grid[1, l], grid[2, l])
      }
    }
  } else {
    Vwm <- array(NA, dim = c(2,2,2))
  }

  Vwm <- ifelse(is.infinite(Vwm), 0, Vwm)

  kappa <- kappa#similparam[7]
  hP0_m0 <- as.vector(modelpriors$hP0_m0)
  hP0_nu0 <- as.double(modelpriors$hP0_nu0)
  hP0_s0 <- as.double(modelpriors$hP0_s0)
  hP0_Lambda0 <- as.double(modelpriors$hP0_Lambda0)

  treatments <- c(as.vector(treatments))

  out <- dm_ppmx_ct_fix(as.integer(iter), as.integer(burn), as.integer(thin),
                 as.integer(nobs), as.vector(treatments), as.integer(PPMx),
                 as.integer(ncon), as.integer(ncat),
                 as.vector(catvec), as.double(kappadp), as.matrix(grid), as.array(Vwm), as.integer(cohesion),
                 as.integer(CC), as.integer(reuse),
                 as.integer(consim), as.integer(similarity),
                 as.integer(calibration), as.integer(coardegree), as.matrix(y), as.matrix(Z),
                 as.matrix(Zpred), as.integer(noprog),
                 as.vector(t(xcon)), as.vector(t(xcat)), as.vector(t(xconp)),
                 as.vector(t(xcatp)), as.integer(npred), as.vector(similparam),
                 as.vector(hP0_m0), as.double(hP0_nu0), as.double(hP0_s0),
                 as.double(hP0_Lambda0), as.integer(update_hierarchy), as.vector(t(beta)),
                 as.integer(hsp), as.vector(mhtunepar), as.integer(A), as.vector(n_a),
                 as.matrix(curr_cluster), as.matrix(card_cluster), as.vector(nclu_curr))

  ###PREPARE OUTPUT
  res <- list()

  #label (clustering)
  label <- list()
  for(t in 1:A){
    label[[t]] <- out$cl_lab[t,,]
  }

  res$label <- label
  res$asstreat <- asstreat
  num_treat <- out$num_treat
  res$num_treat <- num_treat

  #number of cluster
  nclu <- out$nclus
  res$nclu <- nclu#(nT x nout matrix)

  nclu_cs <- apply(nclu, 1, cumsum)#(nt x nout matrix)


  #linear predictor
  eta_out <- out$eta
  res$eta <- eta_out

  #prognostic covariates
  beta <- array(0, dim = c(ncol(Z), ncol(y), nout))

  beta0 <- out$beta
  for(iter in 1:nout){
    for(k in 1:ncol(y)){
      for(q in 1:ncol(Z)){
        h <- q + (k-1) * ncol(Z)
        beta[q, k, iter] <- beta0[h, iter]
      }
    }
  }

  res$beta <- beta
  res$acc_beta <- out$beta_acc
  #pi (multinomial parameter)
  pi_out <- out$pi
  res$pi_out <- pi_out

  res$isypred <- out$yispred
  res$WAIC <- out$WAIC
  res$lpml <- out$lpml
  res$sigmangg <- out$sigma_ngg
  res$kappangg <- out$kappa_ngg

  #posterior predictive
  ypred <- array(0, dim = c(npred, ncol(y), A, nout))
  for(i in 1:nout){
    for(a in 1:A){
      ypred[,,a,i] <- out$ypred[i,1][[1]][,,a]
    }
  }
  res$ypred <- ypred

  pipred <- array(0, dim = c(npred, ncol(y), A, nout))
  for(i in 1:nout){
    for(a in 1:A){
      pipred[,,a,i] <- out$pipred[i,1][[1]][,,a]
    }
  }
  res$pipred <- pipred
  if(any(is.nan(unlist(pipred)))){
    cat("some NaN occurred", "\n")
  }

  res$theta_prior <- out$theta_prior
  res$sigma_prior <- out$sigma_prior
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

