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
  // t - subset of covariates index
  // mm - for auxiliary parameters
  // zi is the latent variable
  // k - categories index
  int l, ll, lll, i, ii, c, p, j, mm, zi, k;
  int dim = y.n_cols;

  int nout = (iter - burn)/(thin); //number of saved iterations

  Rcpp::Rcout << "nout =  " << nout << std::endl;
  //////////////////////////
  //// DM scheme stuff
  //////////////////////////

  arma::vec ypiu(nobs, arma::fill::zeros);
  for(i = 0 ; i < nobs ; i++){
    for(k = 0 ; k < dim ; k++){
      ypiu(i) += y(i, k);
    }
  }

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

  arma::mat mu_star_curr(dim, nobs);
  arma::mat sigma_star_curr(dim * dim, nobs);

  /*for(i = 0; i < nobs; i++){
   mu_star_curr.row(ii) = ran_mvnorm(hP0_m0, hP0_L0, dim).t();
   sigma_star_curr.row(ii) = ran_iwish(hP0_nu0, hP0_V0, dim).t();
  }*/
  mu_star_curr.fill(0); //tbfilled from input
  sigma_star_curr.fill(0);
  int idx = 0;
  for(i = 0; i < dim; i++){
    sigma_star_curr(idx, 0) = 1;
    idx += (dim + 1);
  }
  for(i = 1; i < nobs; i++){
    sigma_star_curr.col(i) = sigma_star_curr.col(0);
  }

  ////////////////////////////////
  //// Model stuff
  ////////////////////////////////

  arma::mat eta(dim, nobs, arma::fill::randu);
  arma::vec eta_flag(nobs, arma::fill::zeros);

  // the linear predictor matrix, represented as a vector, is in the log scale
  // log-linear function on the prognostic marker n x K

  arma::mat loggamma(nobs, dim, arma::fill::zeros);
  for(i = 0 ; i < nobs; i++){
    for(k = 0 ; k < dim; k++){
      loggamma(i, k) = calculate_gamma(eta, curr_clu, k, i, 1);//calculate_gamma(XX, alpha, beta, jj, ii, Log);
    }
  }

  ////////////////////////////////////////////////////
  //// log determinant of matrices stuff
  ////////////////////////////////////////////////////
  int row, col;
  arma::mat work(dim, dim);
  arma::vec wv(dim*dim);
  double ldSig;

  ////////////////////////////////////////////////////
  //// reuse algorithm stuff
  ////////////////////////////////////////////////////

  arma::mat mu_empty(dim, CC, arma::fill::zeros);
  arma::mat sigma_empty(dim * dim, CC, arma::fill::zeros);
  idx = 0;
  for(i = 0; i < dim; i++){
    sigma_empty(idx, 0) = 1;
    idx += (dim + 1);
  }
  for(i = 1; i < CC; i++){
    sigma_empty.col(i) = sigma_empty.col(0);
  }

  /////////////////////////////////////
  // Storage for posterior predictive
  /////////////////////////////////////
  arma::mat ispred_iter(dim, nobs);
  ispred_iter.fill(0.0);

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
  double maxwei, denwei;
  double uu, cweight;
  int newci, id_empty;

  ////////////////////////////////////////
  // Stuff needed parameters
  // (alg 8 Neal 2 step)
  ////////////////////////////////////////

  arma::vec Lnv(dim * dim);
  arma::mat Lnm(dim, dim);
  arma::mat L0mInv(dim, dim);
  arma::mat S0m(dim, dim);
  arma::mat Sigma(dim, dim);
  arma::mat SigmaInv(dim, dim);
  arma::vec mun(dim);
  arma::vec theta(dim);
  arma::vec Snv(dim * dim);
  arma::mat Snm(dim, dim);
  arma::vec etabar(dim);
  etabar.fill(0.0);

  // Stuff to compute lpml (log pseudo marginal likelihood),
  // likelihood, and WAIC widely applicable information criterion (WAIC),
  // also known as Watanabe–Akaike information criterion
  double lpml_iter;
  arma::vec CPOinv(nobs);
  CPOinv.fill(0.0);
  arma::vec like_iter(nobs);
  like_iter.fill(0.0);
  /*double elppdWAIC;
  arma::vec mnlike(nobs);
  mnlike.fill(0.0);
  arma::vec mnllike(nobs);
  mnllike.fill(0.0);*/

  //Stuff for storage and return
  Rcpp::List l_eta_out(2);

  arma::vec nclus(nout);
  arma::mat mu_out(1, dim);
  arma::mat sigma_out(1, dim*dim);
  arma::mat eta_out(1, dim);
  arma::vec Clui(nout * nobs, arma::fill::ones);
  arma::vec like(nout * nobs, arma::fill::ones);
  arma::cube ispred(dim, nobs, nout, arma::fill::zeros);
  ll = 0;
  lll = 0;
  ////////////////////////////////////////
  //
  // HERE COMES THE MCMC
  //
  ////////////////////////////////////////

  for(l = 0; l < iter; l++){
    /*if((l+1) % 10000 == 0){
     Rcpp::Rcout << "mcmc iter =  " << l+1 << std::endl;
     }*/
    //Rcpp::Rcout << "mcmc iter =  " << l+1 << std::endl;
    /* inizializzo latent variables per empty components da P0
     * per il Reuse algorithm (vettori di lunghezza CC)
     */
    if(reuse == 1){
      for(mm = 0; mm < CC; mm++){
        mu_empty.col(mm) = ran_mvnorm(hP0_m0, hP0_L0, dim);
        sigma_empty.col(mm) = ran_iwish(hP0_nu0, hP0_V0, dim);
      }
    }

    for(i = 0; i < nobs; i++){
      /////////////////////////////////////////
      // update the cluster labels with NEAL 8
      /////////////////////////////////////////
      zi = curr_clu(i)-1; //sottraggo 1 perché C conta da 0 ma conviene farlo ogni volta
      //Rcpp::Rcout << "i " << i << ", zi " << zi << std::endl;
      //Rcpp::Rcout << "curr clu adj 0" <<curr_clu.t() << std::endl;
      //Rcpp::Rcout << "nj_curr(zi) " <<nj_curr.t() << std::endl;
      if(nj_curr(zi) > 1){// Case where the nj corresponding to zi is such that nj>1
        nj_curr(zi) -= 1;
      } else {// Case where the nj corresponding to zi is such that nj=1
        dval = R::runif(0.0, 1.0);
        dval *= CC;
        id_empty = floor(dval);//praticamente ho fatto un sample()

        sigma_empty.col(id_empty) = sigma_star_curr.col(zi);
        mu_empty.col(id_empty) = mu_star_curr.col(zi);

        //ADJUST CARDINALITY, \star ptms, NUMB OF CLUSTERS AND LABELS
        nclu_curr -= 1;

        for(ii = 0; ii < nclu_curr; ii++){
          if(ii >= (zi)){
            nj_curr(ii) = nj_curr(ii+1);
            sigma_star_curr.col(ii) = sigma_star_curr.col(ii + 1);
            mu_star_curr.col(ii) = mu_star_curr.col(ii + 1);
          }
        }

        for(ii = 0; ii < nobs; ii++){
          if(curr_clu(ii)>(zi)){
            curr_clu(ii) -= 1;
          }
        }
        //FINE PARTA PROBLEMATICA
      }
      //Rcpp::Rcout << "curr clu adj 1" <<curr_clu.t() << std::endl;
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

          ldSig = logdet(sigma_star_curr.col(j), dim);
          weight(j) = dmvnorm(eta.col(i), mu_star_curr.col(j),
                 sigma_star_curr.col(j), dim, ldSig, 1) + //density
                   log((double) nj_curr(j)) + // cohesion part
                   lgcatY - lgcatN + // Categorical part
                   lgconY - lgconN;  // Continuous part

          if(calibration == 2){
            ldSig = logdet(sigma_star_curr.col(j), dim);
            weight(j) = dmvnorm(eta.col(i), mu_star_curr.col(j),
                   sigma_star_curr.col(j), dim, ldSig, 1)+
                     log((double) nj_curr(j)) + // cohesion part
                     (1/((double)ncon + (double)ncat))*(lgcatY + lgconY - lgcatN - lgconN);
          }
        } else{//quello che segue è PPM (no covariate)
          ldSig = logdet(sigma_star_curr.col(j), dim);
          weight(j) = dmvnorm(eta.col(i), mu_star_curr.col(j),
                 sigma_star_curr.col(j), dim, ldSig, 1)+
                   log((double) nj_curr(j)); // DP cohesion part
        }
      }

      //SIMILARITY & CALIBRATION EMPTY CLUSTERS

      if(reuse == 0){
        for(mm = 0; mm < CC; mm++){
          mu_empty.col(mm) = ran_mvnorm(hP0_m0, hP0_L0, dim);
          sigma_empty.col(mm) = ran_iwish(hP0_nu0, hP0_V0, dim);
        }
      }

      for(j = nclu_curr; j < (nclu_curr + CC); j++){
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

          ldSig = logdet(sigma_empty.col(j - nclu_curr), dim);
          weight(j) = dmvnorm(eta.col(i), mu_empty.col(j - nclu_curr),
                 sigma_empty.col(j - nclu_curr), dim, ldSig, 1) +
                   log(alpha) - log(CC) + //cohesion + auxiliary ptms
                   lgcondraw + // Continuous covariate part
                   lgcatdraw; // categorical covariate part

          if(calibration == 2){
            ldSig = logdet(sigma_empty.col(j - nclu_curr), dim);
            weight(j) = dmvnorm(eta.col(i), mu_empty.col(j - nclu_curr),
                   sigma_empty.col(j - nclu_curr), dim, ldSig, 1) +
                     log(alpha) - log(CC) +
                     (1/((double)ncon + (double)ncat))*(lgcondraw + lgcatdraw);
          }
        } else {// di seguito ppm
          ldSig = logdet(sigma_empty.col(j - nclu_curr), dim);
          weight(j) = dmvnorm(eta.col(i), mu_empty.col(j - nclu_curr),
                 sigma_empty.col(j - nclu_curr), dim, ldSig, 1) +
                   log(alpha) - log(CC); //cohesion + auxiliary ptms
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

          ldSig = logdet(sigma_star_curr.col(j), dim);

          weight(j) = dmvnorm(eta.col(i), mu_star_curr.col(j),
                 sigma_star_curr.col(j), dim, ldSig, 1) +
                   log((double) nj_curr(j)) +  // Cohesion part
                   lgtilYk - lgtilNk; //This takes into account both cont and cat vars
        }

        // calibration for empty clusters
        for(j = nclu_curr; j < nclu_curr + CC; j++){
          ldSig = logdet(sigma_empty.col(j - nclu_curr), dim);
          weight(j) = dmvnorm(eta.col(i), mu_empty.col(j - nclu_curr),
                 sigma_empty.col(j - nclu_curr), dim, ldSig, 1) +
                   log(alpha) - log(CC) +
                   lgtilN(j) - // Continuous covariate part
                   log(sgN);
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
      }else{
        /*id_empty: it identifies which of the augmented has been chosen
         * 2b  page 345 Favaro and Teh, second column
         */
        id_empty = newci - nclu_curr - 1;
        nclu_curr += 1;
        curr_clu(i) = nclu_curr;
        nj_curr(nclu_curr-1) = 1;

        sigma_star_curr.col(nclu_curr-1) = sigma_empty.col(id_empty);
        mu_star_curr.col(nclu_curr-1) = mu_empty.col(id_empty);
        mu_empty.col(id_empty) = ran_mvnorm(hP0_m0, hP0_L0, dim);
        sigma_empty.col(id_empty) = ran_iwish(hP0_nu0, hP0_V0, dim);
      }
    }
    //Rcpp::Rcout << "curr clu fin" <<curr_clu.t() << std::endl;
    //////////////////////////////////////////////////
    // update the cluster value with NEAL 8 (II step)
    //////////////////////////////////////////////////
    for(int row = 0; row < dim; row++){
     for(int col = 0; col < dim; col++){
      L0mInv(row, col) = hP0_L0(row * dim + col);
     }
    }
    L0mInv = arma::inv(L0mInv);

     for(j = 0; j < nclu_curr; j++){
     /*
      * qst non mi serve più perché è come se avessi una sola osservazione.
      * calcolo le full conditional cmq
      *
      *
     arma::mat etatemp(dim, nj_curr(j));
     int idxy = 0;
     for(ii = 0; ii < nobs; ii++){
      if(curr_clu(ii) == (j+1)){
        etatemp.col(idxy) = eta.col(ii);
        etabar += eta.col(ii);
        idxy += 1;
      }
     }

       etabar = etabar/nj_curr(j);

       */

     for(int row = 0; row < dim; row++){
      for(int col = 0; col < dim; col++){
        Sigma(row, col) = sigma_star_curr(row * dim + col, j);
      }
     }
     SigmaInv = arma::inv(Sigma);

            //nSigmaInv = nj_curr(j) * nSigmaInv;
     Lnm = arma::inv(L0mInv + SigmaInv);
     mun = Lnm * (L0mInv * hP0_m0 + SigmaInv);
     for(int row = 0; row < dim; row++){
      for(int col = 0; col < dim; col++){
        Lnv(row * dim + col) = Lnm(row, col);
        S0m(row, col) = hP0_V0(row * dim + col);
        }
      }
     mu_star_curr.col(j) = ran_mvnorm(mun, Lnv, dim);
     /*qui devo calcolarlo attentamente
      *
     Snm = arma::inv(S0m + ((eta.each_col()-mu_star_curr.col(j)) *
       (etatemp.each_col()-mu_star_curr.col(j)).t()));
      */
     Snm = arma::inv(S0m + ((eta.col(j)-mu_star_curr.col(j)) * (eta.col(j)-mu_star_curr.col(j)).t()));
     //Snm = arma::inv(S0m + ((ytemp.each_row()-theta.t()).t() * (ytemp.each_row()-theta.t())));
     for(int row = 0; row < dim; row++){
      for(int col = 0; col < dim; col++){
        Snv(row * dim + col) = Snm(row, col);
      }
     }
     //sigma_star_curr.col(j) = ran_iwish(hP0_nu0+nj_curr(j), Snv, dim);
     sigma_star_curr.col(j) = ran_iwish(hP0_nu0 + 1, Snv, dim);
     }

     for(j = 0; j < nclu_curr; j++){
       l_eta_out = eta_update(JJ, loggamma, nclu_curr, curr_clu, nj_curr,
                  eta.col(j), eta_flag.col(j), mu_star_curr.col(j), sigma_star_curr.col(j), j);
       eta.col(j) = Rcpp::as<arma::vec>(l_eta_out[0]);
       loggamma = Rcpp::as<arma::mat>(l_eta_out[1]);
       eta_flag.col(j) = Rcpp::as<arma::vec>(l_eta_out[2]);
     }

    /*////////////////////////////////////////////////
    * update random variables:
    * - independent gammas for DM sampling scheme JJ, TT
    * - random gamma (trucchetto Raffaele) ss
    ////////////////////////////////////////////////*/

    // update JJ and consequently TT
    for(i = 0 ; i < nobs ; i++){
      TT(ii) = 0.0;
      for(k = 0 ; k < dim; k++){
        JJ(i, k) = R::rgamma(y(i, k)+exp(loggamma(i, k)), 1.0/(ss(i) + 1.0));
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

    /*if((l > (burn-1)) & (l % (thin) == 0)){
     for(i = 0; i < nobs; i++){
     ispred_iter.col(i) = ran_mvnorm(mu_star_curr.col(curr_clu(i)-1),
     sigma_star_curr.col(curr_clu(i)-1), dim);
     }
    }*/

    //////////////////////
    // Save MCMC iterates
    //////////////////////
    if((l > (burn-1)) & ((l + 1) % thin == 0)){
      nclus(ll) = nclu_curr;
      for(i = 0; i < nobs; i++){
        Clui(ll*(nobs) + i) = curr_clu(i);
        like(ll*(nobs) + i) = like_iter(i);
        ispred.slice(ll) = ispred_iter;
      }
      ll += 1;
      for(j = 0; j < nclu_curr; j++){
        mu_out.insert_rows(lll, mu_star_curr.col(j).t());
        sigma_out.insert_rows(lll, sigma_star_curr.col(j).t());
        eta_out.insert_rows(lll, eta.col(j).t());
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
  return Rcpp::List::create(Rcpp::Named("mu") = mu_out,
                            Rcpp::Named("sigma") = sigma_out,
                            Rcpp::Named("cl_lab") = Clui,
                            //Rcpp::Named("like") = like,
                            //Rcpp::Named("ispred") = ispred,
                            Rcpp::Named("nclus") = nclus);//,
                            //Rcpp::Named("WAIC") = WAIC,
                            //Rcpp::Named("lpml") = lpml);
}
