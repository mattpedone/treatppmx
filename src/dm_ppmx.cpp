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
                    arma::vec xcon, arma::vec xcat, arma::vec similparam,
                    arma::vec hP0_m0, arma::vec hP0_L0, double hP0_nu0,
                    arma::vec hP0_V0, arma::vec mhtune){

  // l - MCMC index
  // ll - MCMC index for saving iterates
  // i - individual index
  // ii - second individual index (for double for loops)
  // c - categorical variable index
  // p - number of covariates index
  // j - cluster index
  // jj - empty cluster index
  // t - subset of covariates index
  // mm - for auxiliary parameters
  // zi is the latent variable
  // k - categories index
  int l, ll, lll, i, ii, iii, c, p, j, jj, mm, zi, k;
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

  double max_C, sumx, sumx2, dval;

  max_C = catvec.max(); //maximum number of categories for categorical covariates

  //compute mean and variance of each continuous covariate
  arma::vec xbar(ncon);
  arma::vec s2mle(ncon);

  for(p = 0; p < ncon; p++){
    sumx = 0.0, sumx2 = 0.0;
    for(i = 0; i < nobs; i ++){
      sumx += xcon(i*(ncon) + p);
      sumx2 += xcon(i*(ncon) + p)*xcon(i*(ncon) + p);
    }
    xbar(p) = sumx/((double) nobs);
    s2mle(p) = sumx2/((double) nobs) - xbar(p)*xbar(p);
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
  //// log determinant of matrices stuff
  ////////////////////////////////////////////////////
  //int row, col;
  //arma::mat work(dim, dim);
  //arma::vec wv(dim*dim);
  //double ldSig;

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
  dirweights.fill(similparam(5));
  //	double m0=0.0, s20=0.5, v=1.0, k0=1.0, nu0=1.0;
  double m0 = similparam(0);
  double s20 = similparam(1);
  double v = similparam(2);
  double k0 = similparam(3);
  double nu0 = similparam(4);

  arma::vec xcontmp(nobs);
  arma::vec njc((nobs)*(ncat));

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

  /*arma::vec thetaj(dim, arma::fill::zeros);

  arma::mat Sigma(dim, dim, arma::fill::zeros);
  arma::mat SigmaInv(dim, dim, arma::fill::zeros);

  arma::vec mu0(dim);
  mu0 = hP0_m0;

  arma::vec L0v(dim * dim);
  L0v = hP0_L0;
  arma::mat L0m(dim, dim, arma::fill::zeros);
  arma::mat L0mInv(dim, dim, arma::fill::zeros);

  double nuiw = hP0_nu0;

  arma::vec S0v(dim * dim);
  S0v = hP0_V0;
  arma::mat S0m(dim, dim, arma::fill::zeros);

  arma::vec Sthetav(dim * dim, arma::fill::zeros);
  arma::mat Sthetam(dim, dim, arma::fill::zeros);

  arma::vec Lpv(dim * dim);
  arma::mat Lpm(dim, dim, arma::fill::zeros);
  arma::mat LpmInv(dim, dim, arma::fill::zeros);

  arma::vec mup(dim, arma::fill::zeros);*/

  // Stuff to compute lpml (log pseudo marginal likelihood),
  // likelihood, and WAIC widely applicable information criterion (WAIC),
  // also known as Watanabe–Akaike information criterion
  /*double lpml_iter;
  arma::vec CPOinv(nobs);
  CPOinv.fill(0.0);*/
  arma::vec like_iter(nobs);
  like_iter.fill(0.0);
  /*double elppdWAIC;
  arma::vec mnlike(nobs);
  mnlike.fill(0.0);
  arma::vec mnllike(nobs);
  mnllike.fill(0.0);*/

  ////////////////////////////////////////
  // Stuff for storage and return
  ////////////////////////////////////////

  Rcpp::List l_eta_out(3);

  arma::vec nclus(nout);
  //arma::mat mu_out(1, dim);
  //arma::mat sigma_out(1, dim*dim);
  arma::mat eta_out(1, dim);
  arma::vec Clui(nout * nobs, arma::fill::ones);
  arma::vec like(nout * nobs, arma::fill::ones);
  arma::cube pigreco(nobs, dim, nout, arma::fill::zeros);
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
      } else {// Case where the nj corresponding to zi is such that nj=1
        dval = R::runif(0.0, 1.0);
        dval *= CC;
        id_empty = floor(dval);//praticamente ho fatto un sample()

        eta_empty.col(id_empty) = eta_star_curr.col(zi);

        //ADJUST CARDINALITY, \star ptms, NUMB OF CLUSTERS AND LABELS
        nclu_curr -= 1;

        for(j = 0; j < nclu_curr; j++){
          if(j >= (zi)){
            nj_curr(j) = nj_curr(j+1);
            eta_star_curr.col(j) = eta_star_curr.col(j + 1);
          }
        }

        for(ii = 0; ii < nobs; ii++){
          if(curr_clu(ii)>(zi)){
            curr_clu(ii) -= 1;
          }
        }

        for(ii = 0; ii < nobs; ii++){
          if(ii != i){
            for(k = 0; k < dim; k++){
              loggamma(ii, k) = calculate_gamma(eta_star_curr.col(curr_clu(ii)-1), 0, k, ii, 1);//calculate_gamma(XX, alpha, beta, jj, ii, Log);
            }
          }
        }
        //FINE PARTA PROBLEMATICA
      }

      //SIMILARITY & CALIBRATION CURRENT CLUSTERS
      for(j = 0; j < nclu_curr; j++){
        if(PPMx == 1){
          // Continuous Covariates
          lgconY = 0.0;
          lgconN = 0.0;
          for(p = 0; p < (ncon); p++){
            njtmp = 0;
            sumx = 0.0;
            sumx2 = 0.0;
            for(ii = 0; ii < nobs; ii++){
              if(ii != i){
                if(curr_clu(ii) == j+1){
                  tmp = xcon(ii*(ncon) + p);

                  sumx = sumx + tmp;
                  sumx2 = sumx2 + tmp * tmp;

                  njtmp += 1;
                }
              }
            }
            if(similarity==1){ // Auxilliary
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, sumx, sumx2, xbar(p), njtmp, 0, 0, 1);
                lgconN += lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, xbar(p), s2mle(p), njtmp, 0, 0, 1);
                lgconN += lgcont;
              }
            }
            if(similarity==2){ //Double Dipper
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, sumx, sumx2, xbar(p), njtmp, 1, 0, 1);
                lgconN += lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, xbar(p), s2mle(p), njtmp, 1, 0, 1);
                lgconN += lgcont;
              }
            }


            // now add ith individual back;
            xcontmp(njtmp) = xcon(i*(ncon)+p);
            sumx = sumx + xcon(i*(ncon)+p);
            sumx2 = sumx2 + xcon(i*(ncon)+p)*xcon(i*(ncon)+p);
            njtmp += 1;

            if(similarity==1){ // Auxilliary
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, sumx, sumx2, xbar(p), njtmp, 0, 0, 1);
                lgconY = lgconY + lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, xbar(p), s2mle(p), njtmp, 0, 0, 1);
                lgconY = lgconY + lgcont;
              }
            }
            if(similarity==2){ //Double Dipper
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, sumx, sumx2, xbar(p), njtmp, 1, 0, 1);
                lgconY = lgconY + lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, xbar(p), s2mle(p), njtmp, 1, 0, 1);
                lgconY = lgconY + lgcont;
              }
            }
          }//chiude ciclo su p covariate continue
          // Categorical Covariates
          lgcatY = 0.0;
          lgcatN = 0.0;
          for(p = 0; p < (ncat); p++){
            for(c = 0; c < catvec(p); c++){
              njc(c) = 0;
            }

            njtmp = 0;
            for(ii = 0; ii < nobs; ii++){
              if(ii != i){
                if(curr_clu(ii) == (j + 1)){
                  njc(xcat(ii*(ncat)+p)) = njc(xcat(ii*(ncat)+p)) + 1; // this needs to be a vector
                  njtmp += 1;
                }
              }
            }

            // Auxiliary - Dirichlet-Multinomial
            if(similarity == 1){
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 0, 1);
              lgcatN = lgcatN + lgcatt;
            }
            // Double dipper - Dirichlet-Multinomial
            if(similarity==2){
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 1, 1);
              lgcatN = lgcatN + lgcatt;
            }

            njc(xcat(i*(ncat)+p)) += 1;
            njtmp += 1;

            // Auxiliary - Dirichlet-Multinomial
            if(similarity==1){
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 0, 1);
              lgcatY = lgcatY + lgcatt;
            }
            // Double dipper - Dirichlet-Multinomial
            if(similarity==2){
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 1, 1);
              lgcatY = lgcatY + lgcatt;
            }
          }//chiude ciclo su p covariate discrete

          gtilY(j) = lgconY + lgcatY;
          gtilN(j) = lgconN + lgcatN;

          weight(j) = log((double) nj_curr(j)) + // cohesion part
            lgcatY - lgcatN + // Categorical part
            lgconY - lgconN;  // Continuous part
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
        } else{//quello che segue è PPM (no covariate)
          weight(j) = log((double) nj_curr(j)); // DP cohesion part
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

      for(j = nclu_curr; j < (nclu_curr + CC); j++){
        jj = j - nclu_curr;
        if(PPMx == 1){
          // Continuous Covariates
          lgcondraw = 0.0;
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
          lgcatdraw = 0.0;
          for(p = 0; p < (ncat); p++){
            for(c = 0; c < catvec(p); c++){
              njc(c) = 0;
            }
            njc(xcat(i*(ncat)+p)) = 1;

            // Auxiliary - Dirichlet-Multinomial
            if(similarity == 1){
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 0, 1);
              lgcatdraw += lgcatt;
            }
            // Double dipper - Dirichlet-Multinomial
            if(similarity==2){
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 1, 1);
              lgcatdraw += lgcatt;
            }
          }//chiude ciclo su covariate discrete
          gtilY(j) = lgcondraw + lgcatdraw;
          gtilN(j) = lgcondraw + lgcatdraw;

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

        } else {// di seguito ppm
          weight(j) = log(alpha) - log(CC); //cohesion + auxiliary ptms
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

        sgY=0.0;
        sgN=0.0;
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
          weight(j) = log(alpha) - log(CC) +
            lgtilN(j) - // Continuous covariate part
            log(sgN);
          for(k = 0; k < dim; k++){
            wo = calculate_gamma(eta_empty, jj, k, i, 0);
            weight(j) += wo * log(JJ(i, k)) - lgamma(wo) - JJ(i, k) * (ss(i) + 1.0);;
            if(y(i, k) != 1){
              weight(j) -=  log(JJ(i, k));
            }
          }
        }
      }

      //Rcpp::Rcout << "pesi0: " << weight.t() << std::endl;

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
    }

    //////////////////////////////////////////////////
    // update the cluster value with NEAL 8 (II step)
    //////////////////////////////////////////////////
    /*Sthetam = (eta_star_curr.col(j) - thetaj)*(eta_star_curr.col(j) - thetaj).t();
     Sthetam = arma::inv(S0m + Sthetam);
     Sigma = ran_iwish(nuiw + 1, Sthetam, dim);*/

    for(j = 0; j < nclu_curr; j++){
      l_eta_out = eta_update(JJ, loggamma, nclu_curr, curr_clu, nj_curr,
                eta_star_curr.col(j), eta_flag, hP0_m0, hP0_L0, j);
      eta_star_curr.col(j) = Rcpp::as<arma::vec>(l_eta_out[0]);
      loggamma = Rcpp::as<arma::mat>(l_eta_out[1]);
      eta_flag = Rcpp::as<arma::vec>(l_eta_out[2]);
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

    //////////////////////
    // Save MCMC iterates
    //////////////////////
    if((l > (burn-1)) & ((l + 1) % thin == 0)){
      nclus(ll) = nclu_curr;
      for(i = 0; i < nobs; i++){
        Clui(ll*(nobs) + i) = curr_clu(i);
        like(ll*(nobs) + i) = like_iter(i);
        pigreco.slice(ll).row(i) = JJ.row(i)/TT(i);
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
  // calculate LPML

  /*like_iter(i) = dmvnorm(eta.col(i), mu_star_curr.col(curr_clu(i)-1),
            sigma_star_curr.col(curr_clu(i)-1), dim, ldSig, 0);

  if((l > (burn-1)) & (l % (thin) == 0)){
    // These are needed for WAIC
    //mnlike(i) = mnlike(i) + (like_iter(i))/(double) nout;
    //mnllike(i) = mnllike(i) + log(like_iter(i))/(double) nout;

    CPOinv(i) += pow((double) nout, -1.0) * pow(like_iter(i), -1.0);
  }

  lpml_iter=0.0;

  for(i = 0; i < nobs; i++){
    //lpml_iter += log(pow(CPOinv(i), -1.0));
    lpml_iter += log(pow(CPOinv(i), -1.0));
    //Rcpp::Rcout << "lpml_iter" << lpml_iter <<std::endl;
  }
  double lpml = lpml_iter;*/


  // Computing WAIC  (see Gelman article)

  /*elppdWAIC = 0.0;
   for(i = 0; i < nobs; i++){
   elppdWAIC += (2*mnllike(i) - log(mnlike(i)));
   }

   WAIC = -2*elppdWAIC;*/

  //RETURN
  return Rcpp::List::create(//Rcpp::Named("mu") = mu_out,
                            Rcpp::Named("eta") = eta_out,
                            //Rcpp::Named("sigma") = sigma_out,
                            Rcpp::Named("cl_lab") = Clui,
                            Rcpp::Named("pi") = pigreco,
                            //Rcpp::Named("like") = like,
                            Rcpp::Named("nclus") = nclus);//,
                            //Rcpp::Named("WAIC") = WAIC,
                            //Rcpp::Named("lpml") = lpml);
}
