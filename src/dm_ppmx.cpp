#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <stdlib.h>     /* srand, rand */
#include "utils.h"

// [[Rcpp::export]]
Rcpp::List dm_ppmx(int iter, int burn, int thin, int nobs, int PPMx, int ncon, int ncat,
                   arma::vec catvec, double alpha, int CC, int reuse, int consim,
                   int similarity, int calibration, arma::mat y,
                   arma::vec xcon, arma::vec xcat,
                   arma::vec xconp, arma::vec xcatp, int npred,
                   arma::vec similparam, arma::vec hP0_m0, arma::vec hP0_L0, double hP0_nu0,
                   arma::vec hP0_V0, int upd_hier){

  // l - MCMC index
  // ll - MCMC index for saving iterates
  // i - individual index
  // ii - second individual index (for double for loops)
  // c - categorical variable index
  // p - number of covariates index
  // pp - index for prediction
  // j - cluster index
  // jj - empty cluster index
  // t - subset of covariates index
  // mm - for auxiliary parameters
  // zi is the latent variable
  // k - categories index
  int l, ll, lll, i, ii, iii, c, p, pp, j, jj, mm, zi, k;
  int dim = y.n_cols;

  int nout = (iter - burn)/(thin); //number of saved iterations

  Rcpp::Rcout << "nout =  " << nout << std::endl;
  //////////////////////////
  //// DM scheme stuff
  //////////////////////////

  arma::vec ypiu(nobs, arma::fill::ones);
  /*for(i = 0 ; i < nobs ; i++){
   for(k = 0 ; k < dim ; k++){
   ypiu(i) += y(i, k);
   }
  }*/

  // initialize latent variable: JJ ~ gamma()
  arma::mat JJ(nobs, dim);

  // TT is the vector of normalization constants
  arma::vec TT(nobs, arma::fill::zeros);

  for(i = 0; i < nobs; i++){
    for(k = 0 ; k < dim; k++){
      JJ(i, k) = ((double) y(i, k));
      // a very small number
      if(JJ(i, k) < pow(10.0, -100.0)){
        JJ(i, k) = pow(10.0, -100.0);
      }
      TT(i) += JJ(i, k);
    }
  }

  // initialize the other clever latent variable: ss
  // note that this uses the data to initialize
  arma::vec ss(nobs);
  for(i = 0 ; i < nobs ; i++){
    ss(i) = R::rgamma(ypiu(i), 1.0/TT(i));
  }

  ////////////////////////////////
  //// Predictive covariate stuff
  ////////////////////////////////

  double max_C, sum, sum2, dval;

  max_C = catvec.max(); //maximum number of categories for categorical covariates

  //compute mean and variance of each continuous covariate
  arma::vec xbar(ncon);
  arma::vec s2mle(ncon);

  if(PPMx == 1){
    for(p = 0; p < ncon; p++){
      sum = 0.0, sum2 = 0.0;
      for(i = 0; i < nobs; i++){
        sum += xcon(i*(ncon) + p);
        //Rcpp::Rcout << "idx" << i*(ncon) + p << std::endl;
        sum2 += xcon(i*(ncon) + p)*xcon(i*(ncon) + p);
      }
      //Rcpp::Rcout << "sum" << sum << std::endl;
      xbar(p) = sum/((double) nobs);
      s2mle(p) = sum2/((double) nobs) - xbar(p)*xbar(p);
    }
  }

  //////////////////////////
  ////Cluster-related stuff
  //////////////////////////

  /* inizializzo il vec delle label dei cluster, il vettore delle cardinalità
   * di ciascun cluster e l'int che mi dice il numero di cluster correnti
   */
  arma::vec curr_clu(nobs);
  //each observation is initialized in a separate cluster
  for(i = 0; i < nobs; i++){
    curr_clu(i) = i + 1;//R::rbinom(1, 0.5) + 1;
  }

  arma::vec nj_curr(nobs);
  nj_curr.fill(0);
  //cfr adr1
  for(i = 0; i < nobs; i++){
    for(ii = 0; ii < nobs; ii++){
      if(curr_clu(i) == ii+1) nj_curr(ii) += 1;
    }
  }
  int nclu_curr = 0;
  for(i = 0; i < nobs; i++){
    if(nj_curr(i) > 0) nclu_curr += 1;
  }

  /*
   * questo è il vettore dei pesi, le probabilità per \rho (Alg 8, Neal 2000)
   * Ha lunghezza pari al massimo numeri di componenti possibile + il numero
   * di variabili ausiliarie che alloco permanentemente con il Reuse Alg
   */

  arma::vec weight(nobs + CC);
  arma::vec pweight(nobs + CC);
  double vp = pow(nobs + CC, -1.0);
  weight.fill(vp);

  ////////////////////////////////////////////////////
  //// cluster specific parameters stuff
  ////////////////////////////////////////////////////

  arma::mat eta_star_curr(dim, nobs, arma::fill::randu);

  ////////////////////////////////
  //// Model stuff
  ////////////////////////////////

  arma::vec eta_flag(nobs, arma::fill::zeros);

  // the linear predictor matrix, represented as a vector, is in the log scale
  // log-linear function on the prognostic marker n x K
  int clu_lg;

  arma::mat loggamma(nobs, dim, arma::fill::zeros);
  for(i = 0; i < nobs; i++){
    clu_lg = curr_clu(i)-1;
    for(k = 0; k < dim; k++){
      loggamma(i, k) = calculate_gamma(eta_star_curr, clu_lg, k, i, 1);//calculate_gamma(XX, alpha, beta, jj, ii, Log);
    }
  }

  ////////////////////////////////////////////////////
  //// cluster specific suffiient statistics
  ////////////////////////////////////////////////////

  arma::vec sumx(nobs * ncon, arma::fill::zeros);
  arma::vec sumx2(nobs * ncon, arma::fill::zeros);
  arma::vec njc(nobs * ncat * max_C, arma::fill::zeros);
  arma::vec njctmp(max_C, arma::fill::zeros);
  double auxreal, sumxtmp, sumx2tmp;
  arma::vec auxv(dim);
  int iaux, auxint;

  if(PPMx == 1){
    // Fill in cluster-specific sufficient statistics based on first partition
    for(i = 0; i < nobs; i++){
      for(p = 0; p < ncon; p++){
        sumx((curr_clu(i)-1)*ncon + p) += xcon(i * ncon + p);
        sumx2((curr_clu(i)-1)*ncon + p) += xcon(i * ncon + p) * xcon(i * ncon + p);
      }
      for(p = 0; p < ncat; p++){
        njc(((curr_clu(i)-1)*ncat + p) * max_C + xcat(i * ncat + p)) += 1;
      }
    }
  }

  ////////////////////////////////////////////////////
  //// reuse algorithm stuff
  ////////////////////////////////////////////////////

  arma::mat eta_empty(dim, CC, arma::fill::randu);

  ////////////////////////////////////////
  // Stuff needed for similarities
  ////////////////////////////////////////

  int njtmp;

  double lgconN, lgconY, lgcatN, lgcatY, tmp;
  double lgcont, lgcatt;

  // Similarity function parameters
  // dirichlet denominator parameter
  arma::vec dirweights(max_C);
  dirweights.fill(similparam(6));
  //	double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
  double m0 = similparam(0);
  double s20 = similparam(1);
  double v = similparam(2);
  double k0 = similparam(3);
  double nu0 = similparam(4);
  double n0 = similparam(5);

  //arma::vec xcontmp(nobs);
  //arma::vec njc((nobs)*(ncat));
  //arma::vec njctmp(max_C);

  arma::vec gtilN(nobs + CC);
  gtilN.fill(0.0);
  arma::vec gtilY(nobs + CC);
  gtilY.fill(0.0);
  arma::vec lgtilN(nobs + CC);
  lgtilN.fill(0.0);
  arma::vec lgtilY(nobs + CC);
  lgtilY.fill(0.0);

  double lgcondraw, lgcatdraw;
  double sgY, sgN,  lgtilNk, lgtilYk, maxgtilY, maxgtilN;

  ////////////////////////////////////////
  // Stuf needed for probabilities
  ////////////////////////////////////////

  double wo;
  double maxwei, denwei;
  double uu, cweight;
  int newci, id_empty;

  ////////////////////////////////////////
  // Stuff needed parameters
  // (alg 8 Neal 2 step)
  ////////////////////////////////////////

  arma::vec mu0(dim);
  mu0 = hP0_m0;

  arma::vec L0v(dim * dim);
  L0v = hP0_L0;
  arma::mat L0_mat(dim, dim, arma::fill::zeros);

  double nuiw = hP0_nu0; //scalar parameter for L0-IW

  arma::vec Psi0(dim * dim); //matrix parameter for L0-IW
  Psi0 = hP0_V0;

  arma::vec emme0(dim, arma::fill::zeros); //prior mean for mu0
  arma::vec sigma0(dim * dim, arma::fill::zeros); //prior variance for mu0

  int idx = 0;
  for(i = 0; i < dim; i++){
    sigma0(idx) = 1;
    idx += (dim + 1);
  }

  arma::mat sigma0_mat(dim, dim, arma::fill::zeros);

  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      sigma0_mat(i, j) = sigma0(i*dim + j);
    }
  }

  arma::mat Psi0_mat(dim, dim, arma::fill::zeros);

  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      Psi0_mat(i, j) = Psi0(i*dim + j);
    }
  }

  arma::mat Vwork(dim, dim, arma::fill::zeros); //working matrix for hyperparameters update
  idx = 0;
  for(i = 0; i < dim; i++){
    Vwork(idx) = 1;
    idx += (dim + 1);
  }

  arma::vec Vvecwork(dim*dim);
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      Vvecwork(i*dim + j) = Vwork(i, j);
    }
  }

  arma::mat Swork(dim, dim, arma::fill::zeros); //work
  arma::vec Rwork(dim, arma::fill::zeros); //work
  arma::mat Bwork(dim, dim, arma::fill::zeros);
  arma::vec Bvecwork(dim*dim, arma::fill::zeros);
  arma::mat Awork(dim, dim, arma::fill::zeros);
  arma::vec Avecwork(dim*dim, arma::fill::zeros);

  //Stuff needed for in-sample e out-of-sample prediction
  arma::vec pii(dim);

  // Stuff to compute the posterior predictive
  arma::mat pii_pred(npred, dim, arma::fill::zeros);
  arma::mat ispred(nobs, dim, arma::fill::zeros);

  arma::mat ppred(npred, dim, arma::fill::zeros);
  arma::vec rbpred(npred, arma::fill::zeros);
  arma::vec predclass(npred, arma::fill::zeros);
  arma::vec predclass_prob(npred*(nobs+1), arma::fill::zeros);

  // Stuff to compute lpml (log pseudo marginal likelihood),
  // likelihood, and WAIC widely applicable information criterion (WAIC),
  // also known as Watanabe–Akaike information criterion
  double lpml_iter;
  arma::vec CPOinv(nobs);
  CPOinv.fill(0.0);
  arma::vec like_iter(nobs);
  like_iter.fill(0.0);
  double elppdWAIC, WAIC;
  arma::vec mnlike(nobs);
  mnlike.fill(0.0);
  arma::vec mnllike(nobs);
  mnllike.fill(0.0);

  ////////////////////////////////////////
  // Stuff for storage and return
  ////////////////////////////////////////

  Rcpp::List l_eta_out(3);

  arma::vec nclus(nout);
  //arma::mat mu_out(1, dim);
  //arma::mat sigma_out(1, dim*dim);
  arma::mat eta_out(1, dim);
  arma::vec Clui(nout * nobs, arma::fill::ones);
  arma::vec predclass_out(nout * npred, arma::fill::zeros);
  arma::vec like(nout * nobs, arma::fill::ones);
  arma::cube pigreco(nobs, dim, nout, arma::fill::zeros);
  arma::cube pigreco_pred(npred, dim, nout, arma::fill::zeros);
  arma::cube ispred_out(nobs, dim, nout, arma::fill::ones);
  arma::cube ppred_out(npred, dim, nout, arma::fill::zeros);
  ll = 0;
  lll = 0;
  ////////////////////////////////////////
  //
  // HERE COMES THE MCMC
  //
  ////////////////////////////////////////

  for(l = 0; l < iter; l++){
    /* inizializzo latent variables per empty components da P0
     * per il Reuse algorithm (vettori di lunghezza CC)
     */
    if(reuse == 1){
      for(mm = 0; mm < CC; mm++){
        eta_empty.col(mm) = ran_mvnorm(hP0_m0, hP0_L0, dim);
      }
    }

    for(i = 0; i < nobs; i++){
      /////////////////////////////////////////
      // update the cluster labels with NEAL 8
      /////////////////////////////////////////
      zi = curr_clu(i)-1; //sottraggo 1 perché C conta da 0 ma conviene farlo ogni volta

      if(nj_curr(zi) > 1){// Case where the nj corresponding to zi is such that nj>1
        nj_curr(zi) -= 1;
        if(PPMx == 1){
          // need to reduce the sumx sumx2 to
          for(p = 0; p < ncon; p++){
            sumx(zi * ncon + p) -= xcon(i * ncon + p);
            sumx2(zi * ncon + p) -= xcon(i * ncon + p) * xcon(i * ncon + p);
          }

          // need to reduce the nhc
          for(p = 0; p < ncat; p++){
            njc((zi * ncat + p) * max_C + xcat(i * ncat + p)) -= 1;
          }
        }
      } else {// Case where the nj corresponding to zi is such that nj=1
        //iaux = curr_clu(i);

        dval = R::runif(0.0, 1.0);
        dval *= CC;
        id_empty = floor(dval);//praticamente ho fatto un sample()



        //what follows is Page's relabelling

        iaux = curr_clu(i);
        if(iaux < nclu_curr){
          for(ii = 0; ii < nobs; ii++){
            if(curr_clu(ii) == nclu_curr){
              curr_clu(ii) = iaux;
            }
          }
          curr_clu(i) = nclu_curr;

          auxv = eta_star_curr.col(zi);
          eta_star_curr.col(zi) = eta_empty.col(id_empty);
          eta_empty.col(id_empty) = auxv;

          nj_curr(zi) = nj_curr(nclu_curr-1);
          nj_curr(nclu_curr-1) = 1;

          if(PPMx == 1){
            // need to swap sumx and sumx2
            for(p = 0; p < ncon; p++){
              auxreal = sumx((iaux-1)*(ncon) + p);
              sumx((iaux-1)*(ncon) + p) = sumx((nclu_curr-1)*(ncon) + p);
              sumx((nclu_curr-1)*(ncon) + p) = auxreal;

              auxreal = sumx2((iaux-1)*(ncon) + p);
              sumx2((iaux-1)*(ncon) + p) = sumx2((nclu_curr-1)*(ncon) + p);
              sumx2((nclu_curr-1)*(ncon) + p) = auxreal;
            }

            // need to swap nhc as well
            for(p = 0; p < ncat; p++){
              for(c=0; c<max_C; c++){
                auxint = njc(((iaux-1)*(ncat) + p)*(max_C) + c);
                njc(((iaux-1)*(ncat) + p)*(max_C) + c) = njc(((nclu_curr-1)*(ncat) + p)*(max_C) + c);
                njc(((nclu_curr-1)*(ncat) + p)*(max_C) + c) = auxint;
              }
            }
          }
        }

        nj_curr(nclu_curr-1) -= 1;

        // need to reduce the sumx sumx2
        if(PPMx == 1){
          for(p = 0; p < ncon; p++){
            sumx((nclu_curr-1)*(ncon) + p) -= xcon(i*(ncon)+p);
            sumx2((nclu_curr-1)*(ncon) + p) -= xcon(i*(ncon)+p)*xcon(i*(ncon)+p);
          }

          // need to reduce the nhc
          for(p = 0; p < ncat; p++){
            njc(((nclu_curr-1)*(ncat) + p)*(max_C) + xcat(i*(ncat)+p)) -= 1;
          }
        }

        nclu_curr -= 1;

        /*eta_empty.col(id_empty) = eta_star_curr.col(zi);

        if(PPMx == 1){
          // need to swap sumx and sumx2
          for(p = 0; p < ncon; p++){
            auxreal = sumx(zi*(ncon) + p);
            sumx(zi*(ncon) + p) = sumx((nclu_curr-1)*(ncon) + p);
            sumx((nclu_curr-1)*(ncon) + p) = auxreal;

            auxreal = sumx2(zi*(ncon) + p);
            sumx2(zi*(ncon) + p) = sumx2((nclu_curr-1)*(ncon) + p);
            sumx2((nclu_curr-1)*(ncon) + p) = auxreal;
          }

          // need to swap nhc as well
          for(p = 0; p < ncat; p++){
            for(c=0; c<max_C; c++){
              auxint = njc((zi*(ncat) + p)*(max_C) + c);
              njc((zi*(ncat) + p)*(max_C) + c) = njc(((nclu_curr-1)*(ncat) + p)*(max_C) + c);
              njc(((nclu_curr-1)*(ncat) + p)*(max_C) + c) = auxint;
            }
          }
        }

        //ADJUST CARDINALITY, \star ptms, NUMB OF CLUSTERS AND LABELS
        nclu_curr -= 1;

        for(j = 0; j < nclu_curr; j++){
          if(j >= (zi)){
            nj_curr(j) = nj_curr(j+1);
            eta_star_curr.col(j) = eta_star_curr.col(j + 1);
          }
        }

        if(PPMx == 1){
          for(p = 0; p < ncon; p++){
            sumx((nclu_curr)*(ncon) + p) -= xcon(i*(ncon)+p);
            sumx2((nclu_curr)*(ncon) + p) -= xcon(i*(ncon)+p)*xcon(i*(ncon)+p);
          }

          // need to reduce the nhc
          for(p = 0; p < ncat; p++){
            njc(((nclu_curr)*(ncat) + p)*(max_C) + xcat(i*(ncat)+p)) -= 1;
          }
        }

        for(ii = 0; ii < nobs; ii++){
          if(curr_clu(ii)>(zi)){
            curr_clu(ii) -= 1;
          }
        }*/

        for(ii = 0; ii < nobs; ii++){
          if(ii != i){
            for(k = 0; k < dim; k++){
              loggamma(ii, k) = calculate_gamma(eta_star_curr.col(curr_clu(ii)-1), 0, k, ii, 1);//calculate_gamma(XX, alpha, beta, jj, ii, Log);
            }
          }
        }
        //FINE PARTE PROBLEMATICA
      }

      //SIMILARITY & CALIBRATION CURRENT CLUSTERS
      for(j = 0; j < nclu_curr; j++){
        lgconY = 0.0;
        lgconN = 0.0;
        lgcatY = 0.0;
        lgcatN = 0.0;
        if(PPMx == 1){
          // Continuous Covariates
          for(p = 0; p < (ncon); p++){
            sumxtmp = sumx(j * ncon + p);
            sumx2tmp = sumx2(j * ncon + p);
            //Rcpp::Rcout << "check1: " << nj_curr(j) << std::endl;
            //Rcpp::Rcout << "check2: " << xbar << std::endl;
            if(similarity==1){ // Auxilliary
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(p), nj_curr(j), 0, 0, 1);
                lgconN += lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(p), s2mle(p), nj_curr(j), 0, 0, 1);
                lgconN += lgcont;
              }
            }
            if(similarity==2){ //Double Dipper
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(p), nj_curr(j), 1, 0, 1);
                lgconN += lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(p), s2mle(p), nj_curr(j), 1, 0, 1);
                lgconN += lgcont;
              }
            }

            // now add ith individual back;
            //xcontmp(njtmp) = xcon(i*(ncon)+p);
            sumxtmp += xcon(i*(ncon)+p);
            sumx2tmp += xcon(i*(ncon)+p)*xcon(i*(ncon)+p);
            //njtmp += 1;

            if(similarity==1){ // Auxilliary
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(p), nj_curr(j) + 1, 0, 0, 1);
                lgconY = lgconY + lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(p), s2mle(p), nj_curr(j) + 1, 0, 0, 1);
                lgconY = lgconY + lgcont;
              }
            }
            if(similarity==2){ //Double Dipper
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(p), nj_curr(j) + 1, 1, 0, 1);
                lgconY = lgconY + lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(p), s2mle(p), nj_curr(j) + 1, 1, 0, 1);
                lgconY = lgconY + lgcont;
              }
            }
          }//chiude ciclo su p covariate continue
          // Categorical Covariates
          for(p = 0; p < ncat; p++){
            for(c = 0; c < max_C; c++){
              njctmp(c) = njc((j*ncat+p)*max_C+c);
            }

            // Auxiliary - Dirichlet-Multinomial
            if(similarity == 1){
              lgcatt = gsimcatDM(njctmp, dirweights, catvec(p), 0, 1);
              lgcatN = lgcatN + lgcatt;
            }
            // Double dipper - Dirichlet-Multinomial
            if(similarity==2){
              lgcatt = gsimcatDM(njctmp, dirweights, catvec(p), 1, 1);
              lgcatN = lgcatN + lgcatt;
            }

            njctmp(xcat(i*(ncat)+p)) += 1;
            //njtmp += 1;

            // Auxiliary - Dirichlet-Multinomial
            if(similarity==1){
              lgcatt = gsimcatDM(njctmp, dirweights, catvec(p), 0, 1);
              lgcatY = lgcatY + lgcatt;
            }
            // Double dipper - Dirichlet-Multinomial
            if(similarity==2){
              lgcatt = gsimcatDM(njctmp, dirweights, catvec(p), 1, 1);
              lgcatY = lgcatY + lgcatt;
            }
          }//chiude ciclo su p covariate discrete

          gtilY(j) = lgconY + lgcatY;
          gtilN(j) = lgconN + lgcatN;
          /*Rcpp::Rcout << "lgconY j: " << lgconY << std::endl;
          Rcpp::Rcout << "lgcatY j: " << lgcatY << std::endl;
          Rcpp::Rcout << "lgconN j: " << lgconN << std::endl;
          Rcpp::Rcout << "lgcatN j: " << lgcatN << std::endl;*/
        }// this closes PPMx

        //compute PLAIN cluster probabilities
        weight(j) = log((double) nj_curr(j)) + // cohesion part
          lgcatY - lgcatN + // Categorical part
          lgconY - lgconN;  // Continuous part
        //le similarità le scrivo così invece che con gtilY(j) e gtilN(j), così quando ho PPM, valgono 0 ed è corretto
        //al contrario se è un PPMx lgcatY, lgcatN, lgconY, lgconN hanno i valori del j-esimo cluster
        for(k = 0; k < dim; k++){
          wo = calculate_gamma(eta_star_curr, j, k, i, 0);
          weight(j) += wo * log(JJ(i, k)) - lgamma(wo) - JJ(i, k) * (ss(i) + 1.0);
          if(y(i, k) != 1){
            weight(j) -=  log(JJ(i, k));
          }
        }

        if(calibration == 2){
          weight(j) = log((double) nj_curr(j)) + // cohesion part
            (1/((double)ncon + (double)ncat))*(lgcatY + lgconY - lgcatN - lgconN);
          for(k = 0; k < dim; k++){
            wo = calculate_gamma(eta_star_curr, j, k, i, 0);
            weight(j) += wo * log(JJ(i, k)) - lgamma(wo) - JJ(i, k) * (ss(i) + 1.0);
            if(y(i, k) != 1){
              weight(j) -=  log(JJ(i, k));
            }
          }
        }
      }//chiude la similarity & calibration for current clusters

      //SIMILARITY & CALIBRATION EMPTY CLUSTERS
      if(reuse == 0){
        for(mm = 0; mm < CC; mm++){
          eta_empty.col(mm) = ran_mvnorm(hP0_m0, hP0_L0, dim);
        }
      }

      lgcondraw = 0.0;
      lgcatdraw = 0.0;
      for(j = nclu_curr; j < (nclu_curr + CC); j++){
        jj = j - nclu_curr;
        if(PPMx == 1){
          // Continuous Covariates
          for(p = 0; p < (ncon); p++){
            tmp = xcon(i*(ncon) + p);
            if(similarity==1){ // Auxilliary
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, tmp, tmp*tmp, xbar(p), 1, 0, 0, 1);
                lgcondraw += lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, tmp, tmp*tmp, xbar(p), s2mle(p), 1, 0, 0, 1);
                lgcondraw += lgcont;
              }
            }
            if(similarity==2){ //Double Dipper
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, tmp, tmp*tmp, xbar(p), 1, 1, 0, 1);
                lgcondraw += lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, tmp, tmp*tmp, xbar(p), s2mle(p), 1, 1, 0, 1);
                lgcondraw += lgcont;
              }
            }
          }//chiude ciclo su p covariate continue

          // Categorical Covariates

          for(p = 0; p < (ncat); p++){
            for(c = 0; c < catvec(p); c++){
              njctmp(c) = 0;
            }
            njctmp(xcat(i*(ncat)+p)) = 1;

            // Auxiliary - Dirichlet-Multinomial
            if(similarity == 1){
              lgcatt = gsimcatDM(njctmp, dirweights, catvec(p), 0, 1);
              lgcatdraw += lgcatt;
            }
            // Double dipper - Dirichlet-Multinomial
            if(similarity==2){
              lgcatt = gsimcatDM(njctmp, dirweights, catvec(p), 1, 1);
              lgcatdraw += lgcatt;
            }
          }//chiude ciclo su covariate discrete
          gtilY(j) = lgcondraw + lgcatdraw;
          gtilN(j) = lgcondraw + lgcatdraw;
        }//closes PPMx

        weight(j) = log(alpha) - log(CC) + //cohesion + auxiliary ptms
          lgcondraw + // Continuous covariate part
          lgcatdraw; // categorical covariate part
        for(k = 0; k < dim; k++){
          wo = calculate_gamma(eta_empty, jj, k, i, 0);
          weight(j) += wo * log(JJ(i, k)) - lgamma(wo) - JJ(i, k) * (ss(i) + 1.0);
          if(y(i, k) != 1){
            weight(j) -=  log(JJ(i, k));
          }
        }
        if(calibration == 2){
          weight(j) = log(alpha) - log(CC) +
            (1/((double)ncon + (double)ncat))*(lgcondraw + lgcatdraw);
          for(k = 0; k < dim; k++){
            wo = calculate_gamma(eta_empty, jj, k, i, 0);
            weight(j) += wo * log(JJ(i, k)) - lgamma(wo) - JJ(i, k) * (ss(i) + 1.0);
            if(y(i, k) != 1){
              weight(j) -=  log(JJ(i, k));
            }
          }
        }
      }

      if((calibration == 1) & (PPMx == 1)){
        maxgtilN = gtilN(0);
        maxgtilY = gtilY(0);

        for(j = 1; j < nclu_curr + CC; j++){

          if(maxgtilN < gtilN(j)) maxgtilN = gtilN(j);

          if(j < nclu_curr){
            if(maxgtilY < gtilY(j)) maxgtilY = gtilY(j);
          }
        }

        sgY = 0.0;
        sgN = 0.0;

        for(j = 0; j < nclu_curr + CC; j++){

          lgtilN(j) = gtilN(j) - maxgtilN;
          sgN = sgN + exp(lgtilN(j));

          if(j < nclu_curr){// se è in un cluster non è un signoletto
            lgtilY(j) = gtilY(j) - maxgtilY;
            sgY = sgY + exp(lgtilY(j));
          }
        }
        // Calibrazione prob di cluster esistenti
        for(j = 0; j < nclu_curr; j++){
          lgtilNk = lgtilN(j) - log(sgN);
          lgtilYk = lgtilY(j) - log(sgY);

          weight(j) = log((double) nj_curr(j)) +  // Cohesion part
            lgtilYk - lgtilNk; //cov cont and cat
          for(k = 0; k < dim; k++){
            wo = calculate_gamma(eta_star_curr, j, k, i, 0);
            weight(j) += wo * log(JJ(i, k)) - lgamma(wo) - JJ(i, k) * (ss(i) + 1.0);
            if(y(i, k) != 1){
              weight(j) -=  log(JJ(i, k));
            }
          }
        }

        // calibration for empty clusters
        for(j = nclu_curr; j < nclu_curr + CC; j++){
          jj = j - nclu_curr;
          lgtilNk = lgtilN(j) - log(sgN);
          //lgtilYk = lgtilY(j) - log(sgY);

          weight(j) = log(alpha) - log(CC) +  // Cohesion part
            //lgtilYk -
            lgtilNk;
          //Rcpp::Rcout << "changed" << std::endl;
          /*weight(j) = log(alpha) - log(CC) +
            lgtilN(j) - // Continuous covariate part
            log(sgN);*/
          for(k = 0; k < dim; k++){
            wo = calculate_gamma(eta_empty, jj, k, i, 0);
            weight(j) += wo * log(JJ(i, k)) - lgamma(wo) - JJ(i, k) * (ss(i) + 1.0);
            if(y(i, k) != 1){
              weight(j) -=  log(JJ(i, k));
            }
          }
        }
      }

      //AVOID ZERO IN WEIGHTS
      maxwei = weight(0);
      for(j = 1; j < nclu_curr+ CC; j++){
        if(maxwei < weight(j)) maxwei = weight(j);
      }

      denwei = 0.0;

      for(j = 0; j < nclu_curr + CC; j++){
        weight(j) = exp(weight(j) - maxwei);
        denwei += weight(j);
      }

      for(j = 0; j < nclu_curr + CC; j++){
        pweight(j) = weight(j)/denwei;
        //mysws += pweight(j);
      }

      //sample the new cluster for i-th observation
      uu = R::runif(0.0,1.0);
      cweight = 0.0;
      //newci = id_empty;
      for(j = 0; j < nclu_curr + CC; j++){
        cweight += pweight(j);

        if (uu < cweight){
          newci = j + 1;
          break;
        }
      }
      /*adjust cluster labels and cardinalities
       * in case 1 the i-th subject goes in an existing clsuter
       * in case 2 the i-th subject goes in an empty one (one of the auxiliary)
       */
      if((newci) <= (nclu_curr)){
        //Point 2a page 345 Favaro and Teh, second column
        curr_clu(i) = newci;
        nj_curr(newci - 1) += 1;
        //Rcpp::Rcout << "if 1" << curr_clu.t() << std::endl;

      }else{
        /*id_empty: it identifies which of the augmented has been chosen
         * 2b  page 345 Favaro and Teh, second column
         */
        id_empty = newci - nclu_curr - 1;
        nclu_curr += 1;
        curr_clu(i) = nclu_curr;
        nj_curr(nclu_curr-1) = 1;

        eta_star_curr.col(nclu_curr-1) = eta_empty.col(id_empty);

        //NEED TO UPDATE GAMMA TOO
        for(ii = 0; ii < nobs; ii++){
          for(k = 0; k < dim; k++){
            loggamma(ii, k) = calculate_gamma(eta_star_curr.col(curr_clu(ii)-1), 0, k, ii, 1);//calculate_gamma(XX, alpha, beta, jj, ii, Log);
          }
        }
        eta_empty.col(id_empty) = ran_mvnorm(hP0_m0, hP0_L0, dim);
      }

      if(PPMx == 1){
        // need to now add the xcon to the cluster to which it was assigned;
        for(p = 0; p < ncon; p++){
          sumx((curr_clu(i)-1)*ncon + p) += xcon(i * ncon + p);
          sumx2((curr_clu(i)-1)*ncon + p) += xcon(i * ncon + p) * xcon(i * ncon + p);
        }
        for(p = 0; p < ncat; p++){
          njc(((curr_clu(i)-1)*ncat + p) * max_C + xcat(i * ncat + p)) += 1;
        }
      }
    }

    //////////////////////////////////////////////////
    // update the cluster value with NEAL 8 (II step)
    //////////////////////////////////////////////////
    /*Sthetam = (eta_star_curr.col(j) - thetaj)*(eta_star_curr.col(j) - thetaj).t();
     Sthetam = arma::inv(S0m + Sthetam);
     Sigma = ran_iwish(nuiw + 1, Sthetam, dim);*/

    for(j = 0; j < nclu_curr; j++){
      l_eta_out = eta_update(JJ, loggamma, nclu_curr, curr_clu, nj_curr,
                             eta_star_curr.col(j), eta_flag, mu0, L0v, j);
      eta_star_curr.col(j) = Rcpp::as<arma::vec>(l_eta_out[0]);
      loggamma = Rcpp::as<arma::mat>(l_eta_out[1]);
      eta_flag = Rcpp::as<arma::vec>(l_eta_out[2]);
    }

    if(upd_hier == 1){
      for(j = 0; j < (dim*dim); j++){
        Swork(j) = 0.0;
      }
      //Rcpp::Rcout << "Swork0" << Swork << std::endl;

      for(j = 0; j < nclu_curr; j++){
        Swork += (eta_star_curr.col(j) - emme0) * (eta_star_curr.col(j) - emme0).t();
      }
      //Rcpp::Rcout << "Swork1" << Swork << std::endl;

      Bwork = nuiw*Psi0_mat + Swork;
      Bwork = arma::inv(Bwork);
      for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
          Bvecwork(i*dim + j) = Bwork(i, j);
        }
      }

      //Rcpp::Rcout << "Bwork1" << Bvecwork << std::endl;

      L0v = ran_iwish(nclu_curr + nuiw, Bvecwork, dim);

      for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
          L0_mat(i, j) = L0v(i*dim + j);
        }
      }

      Vwork = arma::inv(sigma0_mat) + nclu_curr*arma::inv(L0_mat);
      Vwork = arma::inv(Vwork);

      for(j = 0; j < dim; j++){
        Rwork(j) = 0;
      }

      for(j = 0; j < nclu_curr; j++){
        Rwork += eta_star_curr.col(j);
      }

      //Rcpp::Rcout << "Rwork" << Rwork.t() << std::endl;

      Awork = Vwork*((arma::inv(sigma0_mat)*emme0)+ arma::inv(L0_mat)*Rwork);
      //Rcpp::Rcout << "here" << std::endl;
      /*for(i = 0; i < dim; i++){
       for(j = 0; j < dim; j++){
       Avecwork(i*dim + j) = Awork(i, j);
       }
      }*/

      for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
          Vvecwork(i*dim + j) = Vwork(i, j);
        }
      }

      mu0 = ran_mvnorm(Awork, Vvecwork, dim);
      //Rcpp::Rcout << "mu0" << mu0.t() << std::endl;
    }

    /*////////////////////////////////////////////////
     * update random variables:
     * - independent gammas for DM sampling scheme JJ, TT
     * - random gamma (trucchetto Raffaele) ss
     ////////////////////////////////////////////////*/

     // update JJ and consequently TT
     for(i = 0; i < nobs; i++){
       TT(i) = 0.0;
       for(k = 0; k < dim; k++){
         if(y(i, k) != 1){
           JJ(i, k) = R::rgamma(y(i, k)+exp(loggamma(i, k)), 1.0/(ss(i) + 1.0));
         } else {
           JJ(i, k) = R::rgamma(y(i, k)+exp(loggamma(i, k) + 1.0), 1.0/(ss(i) + 1.0));
         }
         if(JJ(i, k) < pow(10.0, -100.0)){
           JJ(i, k) = pow(10.0, -100.0);
         }
         TT(i) += JJ(i, k);
       }
     }

     // update latent variables ss
     for(i = 0 ; i < nobs ; i++){
       ss(i) = R::rgamma(ypiu(i), 1.0/TT(i));
     }

     /*////////////////////////////////////////////////
      * in sample prediction to assess model fit
      ////////////////////////////////////////////////*/

      if((l > (burn-1)) & ((l + 1) % thin == 0)){
        for(i = 0; i < nobs; i++){
          for(k = 0; k < dim; k++){
            pii(k) = JJ(i, k)/TT(i);
          }
          ispred.row(i) = rmultinom_rcpp(1, 1, pii);
          like_iter(i) = dmultinom_rcpp(y.row(i).t(), 1, pii, 0);

          //needed for WAIC
          mnlike(i) += like_iter(i)/(double) nout;
          mnllike(i) += log(like_iter(i))/(double) nout;

          //CPO & lpml
          CPOinv(i) += (1/(double) nout)*(1/like_iter(i));
        }
      }
      /*////////////////////////////////////////////////
       * Posterior Predictive
       ////////////////////////////////////////////////*/
       if((l > (burn-1)) & ((l + 1) % thin == 0)){
         for(pp = 0; pp < npred; pp++){
           for(j = 0; j < nclu_curr; j++){

             lgconN=0.0, lgconY=0.0;
             lgcatN=0.0, lgcatY=0.0;

             if(PPMx == 1){
               for(p = 0; p < ncon; p++){
                 sumxtmp = sumx(j * ncon + p);
                 sumx2tmp = sumx2(j * ncon + p);
                 if(similarity==1){ // Auxilliary
                   if(consim==1){//normal normal
                     lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(p), nj_curr(j), 0, 0, 1);
                     lgconN += lgcont;
                   }
                   if(consim==2){//normal normal inverse gamma
                     lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(p), s2mle(p), nj_curr(j), 0, 0, 1);
                     //sufficient statistics da controllare: in Page "sballate"
                     lgconN += lgcont;
                   }
                 }
                 if(similarity==2){ //Double Dipper
                   if(consim==1){//normal normal
                     lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(p), nj_curr(j), 1, 0, 1);
                     //sufficient statistics da controllare: in Page "sballate"
                     lgconN += lgcont;
                   }
                   if(consim==2){//normal normal inverse gamma
                     lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(p), s2mle(p), nj_curr(j), 1, 0, 1);
                     lgconN += lgcont;
                   }
                 }

                 //add the pp-th predictor to cluster
                 sumxtmp += xconp(pp * ncon + p);
                 sumx2tmp += xconp(pp * ncon + p) * xconp(pp * ncon + p);

                 if(similarity==1){ // Auxilliary
                   if(consim==1){//normal normal
                     lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(p), nj_curr(j) + 1, 0, 0, 1);
                     lgconY += lgcont;
                   }
                   if(consim==2){//normal normal inverse gamma
                     lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(p), s2mle(p), nj_curr(j) + 1, 0, 0, 1);
                     //sufficient statistics da controllare: in Page "sballate"
                     lgconY += lgcont;
                   }
                 }
                 if(similarity==2){ //Double Dipper
                   if(consim==1){//normal normal
                     lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(p), nj_curr(j) + 1, 1, 0, 1);
                     //sufficient statistics da controllare: in Page "sballate"
                     lgconY += lgcont;
                   }
                   if(consim==2){//normal normal inverse gamma
                     lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(p), s2mle(p), nj_curr(j) + 1, 1, 0, 1);
                     lgconN += lgcont;
                   }
                 }
               }//this closes the loop on continuous covariates

               for(p = 0; p < ncat; p++){

                 for(c = 0; c < max_C; c++){
                   njctmp(c) = njc[(j*ncat + p)*(max_C) + c];
                 }

                 // Auxiliary - Dirichlet-Multinomial
                 if(similarity == 1){
                   lgcatt = gsimcatDM(njctmp, dirweights, catvec(p), 0, 1);
                   lgcatN = lgcatN + lgcatt;
                 }
                 // Double dipper - Dirichlet-Multinomial
                 if(similarity==2){
                   lgcatt = gsimcatDM(njctmp, dirweights, catvec(p), 1, 1);
                   lgcatN = lgcatN + lgcatt;
                 }

                 njctmp(xcat(pp*(ncat)+p)) += 1;
                 //njtmp += 1;

                 // Auxiliary - Dirichlet-Multinomial
                 if(similarity==1){
                   lgcatt = gsimcatDM(njctmp, dirweights, catvec(p), 0, 1);
                   lgcatY = lgcatY + lgcatt;
                 }
                 // Double dipper - Dirichlet-Multinomial
                 if(similarity==2){
                   lgcatt = gsimcatDM(njctmp, dirweights, catvec(p), 1, 1);
                   lgcatY = lgcatY + lgcatt;
                 }
               }//this closes the loop on categorical covariates

               gtilY(j) = lgconY + lgcatY;
               gtilN(j) = lgconN + lgcatN;

             }//this closes the if on PPMx

             weight(j) = log((double) nj_curr(j)) + // cohesion part
               lgcatY - lgcatN + // Categorical part
               lgconY - lgconN;  // Continuous part

             if(calibration == 2){
               weight(j) = log((double) nj_curr(j)) + // cohesion part
                 (1/((double)ncon + (double)ncat))*(lgcatY + lgconY - lgcatN - lgconN);
             }
           }//this closes the loop on existing clusters

           lgcondraw = 0.0;
           lgcatdraw = 0.0;

           //probabilità che la predittiva metta l'osservazione in un cluster tutto suo
           for(j = nclu_curr; j < (nclu_curr + CC); j++){
             jj = j - nclu_curr;
             if(PPMx == 1){
             // Continuous Covariates
              for(p = 0; p < (ncon); p++){
               tmp = xconp(pp*(ncon) + p);
               if(similarity==1){ // Auxilliary
                 if(consim==1){//normal normal
                   lgcondraw += gsimconNN(m0, v, s20, tmp, tmp*tmp, xbar(p), 1, 0, 0, 1);
                   }
                 if(consim==2){//normal normal inverse gamma
                   lgcondraw += gsimconNNIG(m0, k0, nu0, s20, tmp, tmp*tmp, xbar(p), s2mle(p), 1, 0, 0, 1);
                   }
                 }
               if(similarity==2){ //Double Dipper
                 if(consim==1){//normal normal
                   lgcondraw += gsimconNN(m0, v, s20, tmp, tmp*tmp, xbar(p), 1, 1, 0, 1);
                   }
                 if(consim==2){//normal normal inverse gamma
                   lgcont += gsimconNNIG(m0, k0, nu0, s20, tmp, tmp*tmp, xbar(p), s2mle(p), 1, 1, 0, 1);
                   }
                 }
               }//chiude ciclo su p covariate continue
             // Categorical Covariates

              for(p = 0; p < (ncat); p++){
               for(c = 0; c < catvec(p); c++){
                 njctmp(c) = 0;
               }
               njctmp(xcatp(pp*(ncat)+p)) = 1;

               // Auxiliary - Dirichlet-Multinomial
               if(similarity == 1){
                 lgcatdraw += gsimcatDM(njctmp, dirweights, catvec(p), 0, 1);
               }
               // Double dipper - Dirichlet-Multinomial
               if(similarity==2){
                 lgcatt += gsimcatDM(njctmp, dirweights, catvec(p), 1, 1);
               }
             }//chiude ciclo su covariate discrete
               gtilY(j) = lgcondraw + lgcatdraw;
               gtilN(j) = lgcondraw + lgcatdraw;
              }//closes PPMx

              weight(j) = log(alpha) - log(CC) + //cohesion + auxiliary ptms
               lgcondraw + // Continuous covariate part
               lgcatdraw; // categorical covariate part
              if(calibration == 2){
               weight(j) = log(alpha) - log(CC) +
                 (1/((double)ncon + (double)ncat))*(lgcondraw + lgcatdraw);
              }
          }//chiude loop su empty cluster

           if((calibration == 1) & (PPMx == 1)){
             maxgtilN = gtilN(0);
             maxgtilY = gtilY(0);

             for(j = 1; j < nclu_curr + CC; j++){

               if(maxgtilN < gtilN(j)) maxgtilN = gtilN(j);

               if(j < nclu_curr){
                 if(maxgtilY < gtilY(j)) maxgtilY = gtilY(j);
               }
             }

             sgY = 0.0;
             sgN = 0.0;

             for(j = 0; j < nclu_curr + CC; j++){

               lgtilN(j) = gtilN(j) - maxgtilN;
               sgN = sgN + exp(lgtilN(j));

               if(j < nclu_curr){// If x is included in an existing cluster in cannot be a singleton
                 lgtilY(j) = gtilY(j) - maxgtilY;
                 sgY = sgY + exp(lgtilY(j));
               }
             }
             // Calibrate the unnormalized cluster probabilities
             for(j = 0; j < nclu_curr; j++){
               lgtilNk = lgtilN(j) - log(sgN);
               lgtilYk = lgtilY(j) - log(sgY);

               weight(j) = log((double) nj_curr(j)) +  // Cohesion part
                 lgtilYk - lgtilNk; //This takes into account both cont and cat vars
             }

             // calibration for empty clusters
             for(j = nclu_curr; j < nclu_curr + CC; j++){
               jj = j - nclu_curr;
               lgtilNk = lgtilN(j) - log(sgN);
               //lgtilYk = lgtilY(j) - log(sgY);

               weight(j) = log(alpha) - log(CC) +  // Cohesion part
                 //lgtilYk -
                 lgtilNk;
               //Rcpp::Rcout << "changed" << std::endl;
               /*weight(j) = log(alpha) - log(CC) +
                lgtilN(j) - // Continuous covariate part
                log(sgN);*/
             }
           }//chiude calibrazione 1

           //AVOID ZERO IN WEIGHTS
           maxwei = weight(0);
           for(j = 1; j < nclu_curr+ CC; j++){
             if(maxwei < weight(j)) maxwei = weight(j);
           }

           denwei = 0.0;

           for(j = 0; j < nclu_curr + CC; j++){
             weight(j) = exp(weight(j) - maxwei);
             denwei += weight(j);
           }

           for(j = 0; j < nclu_curr + CC; j++){
             pweight(j) = weight(j)/denwei;
             //mysws += pweight(j);
           }

           //sample the new cluster for i-th observation
           uu = R::runif(0.0,1.0);
           cweight = 0.0;
           //newci = id_empty;
           for(j = 0; j < nclu_curr + CC; j++){
             cweight += pweight(j);

             if (uu < cweight){
               newci = j + 1;
               break;
             }
           }

           /*
            * adjust cluster labels and cardinalities
            */

           arma::vec eta_pred(dim);
           arma::vec loggamma_pred(dim);

           if((newci) <= (nclu_curr)){
             eta_pred = eta_star_curr.col(newci - 1);
           }else{
             eta_pred = ran_mvnorm(mu0, L0v, dim);
           }

           //NEED TO UPDATE GAMMA TOO
           for(k = 0; k < dim; k++){
             loggamma_pred(k) = calculate_gamma(eta_pred, 0, k, ii, 1);
             }

           pii_pred.row(pp) = exp(loggamma_pred).t();

           pii_pred.row(pp) = pii_pred.row(pp)/arma::sum(pii_pred.row(pp));

           ppred.row(pp) = rmultinom_rcpp(1, 1, pii_pred.row(pp).t());
           predclass(pp) = newci;
           //predclass_prob(pp * nobs +j) = weight(j);

         }//this closes the loop on npred covariates
       }//this closes the if for the draws after burnin and thinned

     //////////////////////
     // Save MCMC iterates
     //////////////////////
     if((l > (burn-1)) & ((l + 1) % thin == 0)){
       nclus(ll) = nclu_curr;
       for(i = 0; i < nobs; i++){
         Clui(ll*(nobs) + i) = curr_clu(i);
         like(ll*(nobs) + i) = like_iter(i);
         ispred_out.slice(ll).row(i) = ispred.row(i);
         pigreco.slice(ll).row(i) = JJ.row(i)/TT(i);
       }
       ppred_out.slice(ll) = ppred;
       pigreco_pred.slice(ll) = pii_pred;
       for(pp = 0; pp < npred; pp++){
         predclass_out(ll*npred + pp) = predclass(pp);
       }

       ll += 1;
       for(j = 0; j < nclu_curr; j++){
         //mu_out.insert_rows(lll, mu_star_curr.col(j).t());
         //sigma_out.insert_rows(lll, sigma_star_curr.col(j).t());
         eta_out.insert_rows(lll, eta_star_curr.col(j).t());
         lll += 1;
       }
     }

  }//CLOSES MCMC iterations

   lpml_iter=0.0;

   for(i = 0; i < nobs; i++){
    lpml_iter += log(pow(CPOinv(i), -1.0));
   }
   double lpml = lpml_iter;


  // Computing WAIC  (see Gelman article)

  elppdWAIC = 0.0;
  for(i = 0; i < nobs; i++){
    elppdWAIC += (2*mnllike(i) - log(mnlike(i)));
  }

   WAIC = -2*elppdWAIC;

  //RETURN
  //Rcpp::Rcout << "eta_flag: " << eta_flag.t() << std::endl;
  return Rcpp::List::create(//Rcpp::Named("mu") = mu_out,
    Rcpp::Named("eta") = eta_out,
    Rcpp::Named("eta_acc") = eta_flag,
    Rcpp::Named("cl_lab") = Clui,
    Rcpp::Named("pi") = pigreco,
    Rcpp::Named("pipred") = pigreco_pred,
    Rcpp::Named("yispred") = ispred_out,
    Rcpp::Named("ypred") = ppred_out,
    Rcpp::Named("clupred") = predclass_out,
    //Rcpp::Named("like") = like,
    Rcpp::Named("nclus") = nclus,
    Rcpp::Named("WAIC") = WAIC,
    Rcpp::Named("lpml") = lpml);
}
