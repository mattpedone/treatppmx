#' ppmx prior
#'
#' @param X covariates
#' @param PPMx option for the use of product partition model with covariates. (default is yes)
#' @param cohesion type of cohesion function that is employed for the PPMx prior on partitions. Options are
#'   1 - DirichletProcess-like cohesion (DP) cohesion
#'   2 - Normalized Generalized Gamma Process (NGG) cohesion
#' @param alpha value of \eqn{\alpha} for cohesion function (concentration parameter in DP)
#' @param sigma value of \eqn{\sigma} for cohesion function (reinforcement parameter in NGG)
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
#' @param iter MCMC number of iteration
#' @param burn MCMC iteration discarded due to burnin
#' @param thin thinning for MCMC
#' @param nclu_init number of cluster used for partial correlation initialization for prognostic covariates coefficient (default)
#' @return Number of cluster at each MCMC iteration
#' @export

prior_ppmx <- function(X=NULL, PPMx = 1, cohesion = 2, alpha=1.0, sigma = 0.2,
                          similarity = 1, consim=1, similparam, calibration=0, coardegree = 1,
                          iter=1100, burn=100, thin=1, nclu_init = 5){

  # X - data.frame whose columns are
  # gcontype - similarity function (1 - Auxilliary, 2 - double dipper)
  # consim=1 => that similarity function of continuous covariate is N-N model (v_j is fixed)
  # consim=2 => that similarity functino of continuous covariate is N-NIG model (v_j is unknown)
  # modelPriors = (mu0, s^2, a1, m)
  # simParms = (m0, s20, v2, k0, nu0, dir, alpha)
  # mh - tuning parameters for MCMC updates of sig2 and sig20


  out <- NULL

  if(!is.data.frame(X) & !is.null(X)) X <- data.frame(X)

  cnames <- colnames(X)
  nobs <- nrow(X)

  nout <- (iter-burn)/thin

  if(is.null(X)){
    xcon <- cbind(rep(0,1));
    xcat <- cbind(rep(0,1));
    ncon <- 0
    ncat <- 0
  }

  if(!(is.null(X))){
    nxobs <- ifelse(is.null(X), 0, nrow(X))
    Xall <- X

    # standardize continuous covariates
    if(nxobs > 0){
      #if(sum(!catvars) > 0){
        Xconstd <- Xall[,, drop=FALSE]#apply(Xall[,!catvars, drop=FALSE], 2, scale)#
        xcon <- Xconstd[1:nobs,,drop=FALSE];
        #print(xcon)
        ncon <- ncol(xcon)
    }
    # Now consider the case when number of covariates for prediction are greater than zero
    }

  curr_cluster <- rep(0, nobs)
  card_cluster <- rep(0, nobs)
  work <- kmeans(X, nclu_init, iter.max = 10, nstart = 25)
  for(i in 1:length(work$cluster)){
    curr_cluster[i] <- work$cluster[i]
  }
  for(i in 1:nclu_init){
    card_cluster[i] <- table(work$cluster)[i]
  }
  nclu_curr <- max(work$cluster)

  #curr_cluster <- sample(nclu_init, nobs, T)
  #card_cluster <- c(table(curr_cluster))
  #nclu_curr <- length(card_cluster)
  #calcola 1 colonna V
  #calcola resto matrice
  #passa a funzione cpp
  #   -V
  #   -sigma
  #   -cohesion

  a <- alpha
  Vwm <- matrix(0, nrow = nobs+1, ncol = nobs+1)
  if(cohesion == 2){
    sigma <- sigma
    Vwm <- matrix(NA, nrow = nobs+1, ncol = nobs+1)
    Vwm[1, 1] <- 1

    for(n in 2:(nobs+1)){
      Vwm[n, (1:n)] <- vweights::computev(n, sigma, a)
    }

    Vwm <- log(Vwm)
  }

  alpha <- alpha
  CC <- 1
  ncat <- 0

  out <- prior_ppmx_core(as.integer(iter), as.integer(burn), as.integer(thin),
                    as.integer(nobs), as.integer(PPMx),
                    as.integer(ncon), as.integer(ncat), as.double(alpha), as.double(sigma), as.matrix(Vwm),
                    as.integer(cohesion),
                    as.integer(CC),
                    as.integer(consim), as.integer(similarity),
                    as.integer(calibration), as.integer(coardegree),
                    as.vector(t(xcon)), as.vector(similparam),
                    as.vector(curr_cluster), as.vector(card_cluster), as.integer(nclu_curr))

  ###PREPARE OUTPUT
  res <- setNames(tabulate(out, nbins = nobs), 1:nobs)

  return <- res

}
