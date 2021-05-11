#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <stdlib.h>     /* srand, rand */
#include "utils_ct.h"

// [[Rcpp::export]]
Rcpp::List dm_ppmx(int iter, int burn, int thin, int nobs, arma::vec treatments, int PPMx, int ncon, int ncat,
                   arma::vec catvec, double alpha, int CC, int reuse, int consim,
                   int similarity, int calibration, arma::mat y, arma::mat z,
                   arma::mat zpred, arma::vec xcon, arma::vec xcat,
                   arma::vec xconp, arma::vec xcatp, int npred,
                   arma::vec similparam, arma::vec hP0_m0, arma::vec hP0_L0, double hP0_nu0,
                   arma::vec hP0_V0, int upd_hier, arma::vec initbeta, int hsp,
                   arma::vec mhtunepar){

  // l - MCMC index
  // ll - MCMC index for saving iterates
  // i - individual index
  // ii - second individual index (for double for loops)
  // c - categorical variable index
  // p - number of covariates index
  // pp - index for subjects in the prediction set
  // j - cluster index
  // jj - empty cluster index
  // mm - for auxiliary parameters
  // zi is the latent variable
  // k - categories index
  // q - prognostico covariate index
  // h - prognostico covariate "instrumental" index
  // tt - treatments index
  int l, ll, lll, i, ii, iii, c, p, pp, j, jj, mm, zi, k, q, h, tt;
  int dim = y.n_cols;
  int Q = z.n_cols;

  int nout = (iter - burn)/(thin); //number of saved iterations

  //Rcpp::Rcout << "nout =  " << nout << std::endl;
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

  //////////////////////////
  ////Multiple Treatments stuff
  //////////////////////////

  int nT = treatments.max(); //numero di trattamenti
  arma::vec num_treat(nT, arma::fill::zeros);// numerosità x trattamento (table(tratments))

  //tt = 0;
  for(i = 0; i < nobs; i++){
    for(tt = 0; tt < nT; tt++){
      if(treatments(i) == tt){
        num_treat(tt) += 1;
      }
    }
  }

  ////////////////////////////////
  //// mean and variance of each continuous covariate
  ////////////////////////////////

  double max_C, dval;
  max_C = catvec.max(); //maximum number of categories for categorical covariates

  arma::vec sum(nT, arma::fill::zeros);
  arma::vec sum2(nT, arma::fill::zeros);

  arma::mat xbar(nT, ncon);
  arma::mat s2mle(nT, ncon);

  if(PPMx == 1){
    for(tt = 0; tt < nT; tt++){
      for(p = 0; p < ncon; p++){
        for(i = 0; i < nobs; i++){
          if(treatments(i) == tt){
            sum(tt) += xcon(i*(ncon) + p);
            sum2(tt) += xcon(i*(ncon) + p)*xcon(i*(ncon) + p);
          }
        }
        //Rcpp::Rcout << "sum" << sum << std::endl;
        xbar(tt, p) = sum(tt)/((double) num_treat(tt));
        s2mle(tt, p) = sum2(tt)/((double) num_treat(tt)) - xbar(tt, p)*xbar(tt, p);
      }
    }
  }

  //////////////////////////
  ////Cluster-related stuff
  //////////////////////////

  /* inizializzo il vec delle label dei cluster, il vettore delle cardinalità
   * di ciascun cluster e l'int che mi dice il numero di cluster correnti
   */
  arma::mat curr_clu(nT, num_treat.max());
  //each observation is initialized in a separate cluster
  for(tt = 0; tt < nT; tt++){
    for(i = 0; i < num_treat(tt); i++){
      curr_clu(tt, i) = i + 1;//R::rbinom(1, 0.5) + 1;
    }
  }

  arma::mat nj_curr(nT, num_treat.max());//numerosità in ciascun cluster x ciascun trattamento
  nj_curr.fill(0);
  //cfr adr1
  for(tt = 0; tt < nT; tt++){
    for(i = 0; i < num_treat(tt); i++){
      for(ii = 0; ii < num_treat(tt); ii++){
        if(curr_clu(tt, i) == ii+1) nj_curr(tt, ii) += 1;
      }
    }
  }
  arma::vec nclu_curr(nT, arma::fill::zeros);//numero di cluster in ciascun trattamento
  for(tt = 0; tt < nT; tt++){
    for(i = 0; i < num_treat(tt); i++){
      if(nj_curr(tt, i) > 0) nclu_curr(tt) += 1;
    }
  }

  /*
   * questo è il vettore dei pesi, le probabilità per \rho (Alg 8, Neal 2000)
   * Ha lunghezza pari al massimo numeri di componenti possibile + il numero
   * di variabili ausiliarie che alloco permanentemente con il Reuse Alg
   */

  arma::mat weight(nT, num_treat.max() + CC);
  arma::mat pweight(nT, num_treat.max() + CC);
  double vp;
  for(tt = 0; tt < nT; tt++){
    vp = pow(num_treat(tt) + CC, -1.0);
    weight.row(tt).fill(vp);
  }
  //controlla se weight è riempito bene

  ////////////////////////////////////////////////////
  //// cluster specific parameters stuff
  ////////////////////////////////////////////////////

  arma::cube eta_star_curr(dim, num_treat.max(), nT, arma::fill::randu);

  ////////////////////////////////
  //// Model stuff
  ////////////////////////////////

  arma::mat eta_flag(nT, num_treat.max(), arma::fill::zeros);
  arma::vec sumtotclu(nT, arma::fill::zeros);

  arma::vec beta = initbeta;
  arma::vec beta_temp(Q, arma::fill::zeros);
  arma::mat beta_flag(Q, dim, arma::fill::zeros);
  double mu_beta = 0.0;
  arma::vec sigma_beta(Q*dim, arma::fill::ones);
  arma::vec lambda(Q*dim, arma::fill::ones);
  double tau = 1.0;

  // the linear predictor matrix, is in the log scale
  // log-linear function on the prognostic marker n x K
  int clu_lg;

  arma::mat loggamma(nobs, dim, arma::fill::zeros);
  for(tt = 0; tt < nT; tt++){
    for(i = 0; i < nobs; i++){
      if(treatments(i) == tt){
        clu_lg = curr_clu(tt, i)-1;
        for(k = 0; k < dim; k++){
          loggamma(i, k) = calculate_gamma(eta_star_curr.slice(tt), z, beta, clu_lg, k, i, 1);
        }
      }
    }
  }

  double opz1, opz2;
  ////////////////////////////////////////////////////
  //// cluster specific suffiient statistics
  ////////////////////////////////////////////////////

  arma::mat sumx(nT, nobs * ncon, arma::fill::zeros);
  arma::mat sumx2(nT, nobs * ncon, arma::fill::zeros);
  arma::mat njc(nT, nobs * ncat * max_C, arma::fill::zeros);
  arma::vec njctmp(max_C, arma::fill::zeros);
  double auxreal, sumxtmp, sumx2tmp;

  arma::vec auxv(dim);
  int iaux, auxint;

  if(PPMx == 1){
    // Fill in cluster-specific sufficient statistics based on first partition
    for(tt = 0; tt < nT; tt++){
      for(i = 0; i < nobs; i++){
        if(treatments(i) == tt){
          for(p = 0; p < ncon; p++){
            sumx(tt, (curr_clu(i)-1)*ncon + p) += xcon(i * ncon + p);
            sumx2(tt, (curr_clu(i)-1)*ncon + p) += xcon(i * ncon + p) * xcon(i * ncon + p);
          }
          for(p = 0; p < ncat; p++){
            njc(tt, ((curr_clu(i)-1)*ncat + p) * max_C + xcat(i * ncat + p)) += 1;
          }
        }
      }
    }
  }

  ////////////////////////////////////////////////////
  //// reuse algorithm stuff
  ////////////////////////////////////////////////////

  arma::cube eta_empty(dim, CC, nT, arma::fill::randu);

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

  arma::vec gtilN(num_treat.max() + CC);
  gtilN.fill(0.0);
  arma::vec gtilY(num_treat.max() + CC);
  gtilY.fill(0.0);
  arma::vec lgtilN(num_treat.max() + CC);
  lgtilN.fill(0.0);
  arma::vec lgtilY(num_treat.max() + CC);
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
  arma::mat eta_pred(nT, dim);
  arma::vec loggamma_pred(dim);

  arma::cube pii_pred(npred, dim, nT, arma::fill::zeros);
  arma::mat ispred(nobs, dim, arma::fill::zeros);

  arma::cube ppred(npred, dim, nT, arma::fill::zeros);
  arma::vec rbpred(npred, arma::fill::zeros);
  arma::mat predclass(npred, nT, arma::fill::zeros);
  //arma::vec predclass_prob(npred*(nobs+1), arma::fill::zeros);

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
  Rcpp::List l_beta_out(3);

  arma::vec nclus(nout);
  //arma::mat sigma_out(1, dim*dim);
  arma::mat eta_out(1, dim);
  arma::mat beta_out(Q * dim, nout, arma::fill::ones);
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
    for(tt = 0 ; tt < nT; tt++){
      /* inizializzo latent variables per empty components da P0
       * per il Reuse algorithm (vettori di lunghezza CC)
       */
      if(reuse == 1){
        for(mm = 0; mm < CC; mm++){
          eta_empty.slice(tt).col(mm) = ran_mvnorm(hP0_m0, hP0_L0, dim);
        }
      }

      for(i = 0; i < nobs; i++){
        if(treatments(i) == tt){
          /////////////////////////////////////////
          // update the cluster labels with NEAL 8
          /////////////////////////////////////////
          zi = curr_clu(tt, i)-1; //sottraggo 1 perché C conta da 0 ma conviene farlo ogni volta

          if(nj_curr(tt, zi) > 1){// Case where the nj corresponding to zi is such that nj>1
            nj_curr(tt, zi) -= 1;
            if(PPMx == 1){
              // need to reduce the sumx sumx2 to
              for(p = 0; p < ncon; p++){
                sumx(tt, zi * ncon + p) -= xcon(i * ncon + p);
                sumx2(tt, zi * ncon + p) -= xcon(i * ncon + p) * xcon(i * ncon + p);
              }

              // need to reduce the nhc
              for(p = 0; p < ncat; p++){
                njc(tt, (zi * ncat + p) * max_C + xcat(i * ncat + p)) -= 1;
              }
            }
          } else {// Case where the nj corresponding to zi is such that nj=1
            //iaux = curr_clu(i);

            dval = R::runif(0.0, 1.0);
            dval *= CC;
            id_empty = floor(dval);//praticamente ho fatto un sample()

            //what follows is Page's relabelling

            iaux = curr_clu(tt, i);
            if(iaux < nclu_curr(tt)){
              for(ii = 0; ii < num_treat(tt); ii++){
                if(curr_clu(tt, ii) == nclu_curr(tt)){
                  curr_clu(tt, ii) = iaux;
                }
              }
              curr_clu(tt, i) = nclu_curr(tt);

              auxv = eta_star_curr.slice(tt).col(zi);
              eta_star_curr.slice(tt).col(zi) = eta_empty.slice(tt).col(id_empty);
              eta_empty.slice(tt).col(id_empty) = auxv;

              nj_curr(tt, zi) = nj_curr(tt, nclu_curr(tt)-1);
              nj_curr(tt, nclu_curr(tt)-1) = 1;

              if(PPMx == 1){
                // need to swap sumx and sumx2
                for(p = 0; p < ncon; p++){
                  auxreal = sumx(tt, (iaux-1)*(ncon) + p);
                  sumx(tt, (iaux-1)*ncon + p) = sumx(tt, (nclu_curr(tt)-1)*(ncon) + p);
                  sumx(tt, (nclu_curr(tt)-1)*(ncon) + p) = auxreal;

                  auxreal = sumx2(tt, (iaux-1)*(ncon) + p);
                  sumx2(tt, (iaux-1)*(ncon) + p) = sumx2(tt, (nclu_curr(tt)-1)*(ncon) + p);
                  sumx2(tt, (nclu_curr(tt)-1)*(ncon) + p) = auxreal;
                }

                // need to swap nhc as well
                for(p = 0; p < ncat; p++){
                  for(c=0; c<max_C; c++){
                    auxint = njc(tt, ((iaux-1)*(ncat) + p)*(max_C) + c);
                    njc(tt, ((iaux-1)*(ncat) + p)*(max_C) + c) = njc(tt, ((nclu_curr(tt)-1)*(ncat) + p)*(max_C) + c);
                    njc(tt, ((nclu_curr(tt)-1)*(ncat) + p)*(max_C) + c) = auxint;
                  }
                }
              }
            }

            nj_curr(tt, nclu_curr(tt)-1) -= 1;

            // need to reduce the sumx sumx2
            if(PPMx == 1){
              for(p = 0; p < ncon; p++){
                sumx(tt, (nclu_curr(tt)-1)*(ncon) + p) -= xcon(i*(ncon)+p);
                sumx2(tt, (nclu_curr(tt)-1)*(ncon) + p) -= xcon(i*(ncon)+p)*xcon(i*(ncon)+p);
              }

              // need to reduce the nhc
              for(p = 0; p < ncat; p++){
                njc(tt, ((nclu_curr(tt)-1)*(ncat) + p)*(max_C) + xcat(i*(ncat)+p)) -= 1;
              }
            }

            nclu_curr(tt) -= 1;

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
              if(treatments(ii) == tt){
                if(ii != i){
                  clu_lg = curr_clu(tt, ii) -1;
                  for(k = 0; k < dim; k++){
                    loggamma(ii, k) = calculate_gamma(eta_star_curr.slice(tt).col(clu_lg),
                             z, beta, 0, k, ii, 1);
                  }
                }
              }
            }
            //FINE PARTE PROBLEMATICA
          }

          //SIMILARITY & CALIBRATION CURRENT CLUSTERS
          for(j = 0; j < nclu_curr(tt); j++){
            lgconY = 0.0;
            lgconN = 0.0;
            lgcatY = 0.0;
            lgcatN = 0.0;
            if(PPMx == 1){
              // Continuous Covariates
              for(p = 0; p < (ncon); p++){
                sumxtmp = sumx(tt, j * ncon + p);
                sumx2tmp = sumx2(tt, j * ncon + p);
                //Rcpp::Rcout << "check1: " << nj_curr(j) << std::endl;
                //Rcpp::Rcout << "check2: " << xbar << std::endl;
                if(similarity==1){ // Auxilliary
                  if(consim==1){//normal normal
                    lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(tt, p), nj_curr(tt, j), 0, 0, 1);
                    lgconN += lgcont;
                  }
                  if(consim==2){//normal normal inverse gamma
                    lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(tt, p), s2mle(tt, p), nj_curr(tt, j), 0, 0, 1);
                    lgconN += lgcont;
                  }
                }
                if(similarity==2){ //Double Dipper
                  if(consim==1){//normal normal
                    lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(tt, p), nj_curr(tt, j), 1, 0, 1);
                    lgconN += lgcont;
                  }
                  if(consim==2){//normal normal inverse gamma
                    lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(tt, p), s2mle(tt, p), nj_curr(tt, j), 1, 0, 1);
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
                    lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(tt, p), nj_curr(tt, j) + 1, 0, 0, 1);
                    lgconY = lgconY + lgcont;
                  }
                  if(consim==2){//normal normal inverse gamma
                    lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(tt, p), s2mle(tt, p), nj_curr(tt, j) + 1, 0, 0, 1);
                    lgconY = lgconY + lgcont;
                  }
                }
                if(similarity==2){ //Double Dipper
                  if(consim==1){//normal normal
                    lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(tt, p), nj_curr(tt, j) + 1, 1, 0, 1);
                    lgconY = lgconY + lgcont;
                  }
                  if(consim==2){//normal normal inverse gamma
                    lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(tt, p), s2mle(tt, p), nj_curr(tt, j) + 1, 1, 0, 1);
                    lgconY = lgconY + lgcont;
                  }
                }
              }//chiude ciclo su p covariate continue
              // Categorical Covariates
              for(p = 0; p < ncat; p++){
                for(c = 0; c < max_C; c++){
                  njctmp(c) = njc(tt, (j*ncat+p)*max_C+c);
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
            weight(tt, j) = log((double) nj_curr(tt, j)) + // cohesion part
              lgcatY - lgcatN + // Categorical part
              lgconY - lgconN;  // Continuous part
            //le similarità le scrivo così invece che con gtilY(j) e gtilN(j), così quando ho PPM, valgono 0 ed è corretto
            //al contrario se è un PPMx lgcatY, lgcatN, lgconY, lgconN hanno i valori del j-esimo cluster
            for(k = 0; k < dim; k++){
              wo = calculate_gamma(eta_star_curr.slice(tt), z, beta, j, k, i, 0);
              weight(tt, j) += wo * log(JJ(i, k)) - lgamma(wo) - JJ(i, k) * (ss(i) + 1.0);
              if(y(i, k) != 1){
                weight(tt, j) -=  log(JJ(i, k));
              }
            }

            if(calibration == 2){
              weight(tt, j) = log((double) nj_curr(tt, j)) + // cohesion part
                (1/((double)ncon + (double)ncat))*(lgcatY + lgconY - lgcatN - lgconN);
              for(k = 0; k < dim; k++){
                wo = calculate_gamma(eta_star_curr.slice(tt), z, beta, j, k, i, 0);
                weight(tt, j) += wo * log(JJ(i, k)) - lgamma(wo) - JJ(i, k) * (ss(i) + 1.0);
                if(y(i, k) != 1){
                  weight(tt, j) -=  log(JJ(i, k));
                }
              }
            }
          }//chiude la similarity & calibration for current clusters

          //SIMILARITY & CALIBRATION EMPTY CLUSTERS
          if(reuse == 0){
            for(mm = 0; mm < CC; mm++){
              eta_empty.slice(tt).col(mm) = ran_mvnorm(hP0_m0, hP0_L0, dim);
            }
          }

          lgcondraw = 0.0;
          lgcatdraw = 0.0;
          for(j = nclu_curr(tt); j < (nclu_curr(tt) + CC); j++){
            jj = j - nclu_curr(tt);
            if(PPMx == 1){
              // Continuous Covariates
              for(p = 0; p < (ncon); p++){
                tmp = xcon(i*(ncon) + p);
                if(similarity==1){ // Auxilliary
                  if(consim==1){//normal normal
                    lgcont = gsimconNN(m0, v, s20, tmp, tmp*tmp, xbar(tt, p), 1, 0, 0, 1);
                    lgcondraw += lgcont;
                  }
                  if(consim==2){//normal normal inverse gamma
                    lgcont = gsimconNNIG(m0, k0, nu0, s20, tmp, tmp*tmp, xbar(tt, p), s2mle(tt, p), 1, 0, 0, 1);
                    lgcondraw += lgcont;
                  }
                }
                if(similarity==2){ //Double Dipper
                  if(consim==1){//normal normal
                    lgcont = gsimconNN(m0, v, s20, tmp, tmp*tmp, xbar(tt, p), 1, 1, 0, 1);
                    lgcondraw += lgcont;
                  }
                  if(consim==2){//normal normal inverse gamma
                    lgcont = gsimconNNIG(m0, k0, nu0, s20, tmp, tmp*tmp, xbar(tt, p), s2mle(tt, p), 1, 1, 0, 1);
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

            weight(tt, j) = log(alpha) - log(CC) + //cohesion + auxiliary ptms
              lgcondraw + // Continuous covariate part
              lgcatdraw; // categorical covariate part
            for(k = 0; k < dim; k++){
              wo = calculate_gamma(eta_empty.slice(tt), z, beta, jj, k, i, 0);
              weight(tt, j) += wo * log(JJ(i, k)) - lgamma(wo) - JJ(i, k) * (ss(i) + 1.0);
              if(y(i, k) != 1){
                weight(tt, j) -=  log(JJ(i, k));
              }
            }
            if(calibration == 2){
              weight(j) = log(alpha) - log(CC) +
                (1/((double)ncon + (double)ncat))*(lgcondraw + lgcatdraw);
              for(k = 0; k < dim; k++){
                wo = calculate_gamma(eta_empty.slice(tt), z, beta, jj, k, i, 0);
                weight(tt, j) += wo * log(JJ(i, k)) - lgamma(wo) - JJ(i, k) * (ss(i) + 1.0);
                if(y(i, k) != 1){
                  weight(tt, j) -=  log(JJ(i, k));
                }
              }
            }
          }

          if((calibration == 1) & (PPMx == 1)){
            maxgtilN = gtilN(0);
            maxgtilY = gtilY(0);

            for(j = 1; j < nclu_curr(tt) + CC; j++){

              if(maxgtilN < gtilN(j)) maxgtilN = gtilN(j);

              if(j < nclu_curr(tt)){
                if(maxgtilY < gtilY(j)) maxgtilY = gtilY(j);
              }
            }

            sgY = 0.0;
            sgN = 0.0;

            for(j = 0; j < nclu_curr(tt) + CC; j++){

              lgtilN(j) = gtilN(j) - maxgtilN;
              sgN = sgN + exp(lgtilN(j));

              if(j < nclu_curr(tt)){// se è in un cluster non è un signoletto
                lgtilY(j) = gtilY(j) - maxgtilY;
                sgY = sgY + exp(lgtilY(j));
              }
            }
            // Calibrazione prob di cluster esistenti
            for(j = 0; j < nclu_curr(tt); j++){
              lgtilNk = lgtilN(j) - log(sgN);
              lgtilYk = lgtilY(j) - log(sgY);

              weight(tt, j) = log((double) nj_curr(tt, j)) +  // Cohesion part
                lgtilYk - lgtilNk; //cov cont and cat
              for(k = 0; k < dim; k++){
                wo = calculate_gamma(eta_star_curr.slice(tt), z, beta, j, k, i, 0);
                weight(tt, j) += wo * log(JJ(i, k)) - lgamma(wo) - JJ(i, k) * (ss(i) + 1.0);
                if(y(i, k) != 1){
                  weight(tt, j) -=  log(JJ(i, k));
                }
              }
            }

            // calibration for empty clusters
            for(j = nclu_curr(tt); j < nclu_curr(tt) + CC; j++){
              jj = j - nclu_curr(tt);
              lgtilNk = lgtilN(j) - log(sgN);
              //lgtilYk = lgtilY(j) - log(sgY);

              weight(tt, j) = log(alpha) - log(CC) +  // Cohesion part
                //lgtilYk -
                lgtilNk;
              //Rcpp::Rcout << "changed" << std::endl;
              /*weight(j) = log(alpha) - log(CC) +
               lgtilN(j) - // Continuous covariate part
               log(sgN);*/
              for(k = 0; k < dim; k++){
                wo = calculate_gamma(eta_empty.slice(tt), z, beta, jj, k, i, 0);
                weight(tt, j) += wo * log(JJ(i, k)) - lgamma(wo) - JJ(i, k) * (ss(i) + 1.0);
                if(y(i, k) != 1){
                  weight(tt, j) -=  log(JJ(i, k));
                }
              }
            }
          }

          //AVOID ZERO IN WEIGHTS
          maxwei = weight(tt, 0);
          for(j = 1; j < nclu_curr(tt)+ CC; j++){
            if(maxwei < weight(tt, j)) maxwei = weight(tt, j);
          }

          denwei = 0.0;

          for(j = 0; j < nclu_curr(tt) + CC; j++){
            weight(tt, j) = exp(weight(tt, j) - maxwei);
            denwei += weight(tt, j);
          }

          for(j = 0; j < nclu_curr(tt) + CC; j++){
            pweight(tt, j) = weight(tt, j)/denwei;
            //mysws += pweight(j);
          }

          //sample the new cluster for i-th observation
          uu = R::runif(0.0,1.0);
          cweight = 0.0;
          //newci = id_empty;
          for(j = 0; j < nclu_curr(tt) + CC; j++){
            cweight += pweight(tt, j);

            if (uu < cweight){
              newci = j + 1;
              break;
            }
          }
          /*adjust cluster labels and cardinalities
           * in case 1 the i-th subject goes in an existing clsuter
           * in case 2 the i-th subject goes in an empty one (one of the auxiliary)
           */
          if((newci) <= (nclu_curr(tt))){
            //Point 2a page 345 Favaro and Teh, second column
            curr_clu(tt, i) = newci;
            nj_curr(tt, newci - 1) += 1;
            //Rcpp::Rcout << "if 1" << curr_clu.t() << std::endl;
            for(k = 0; k < dim; k++){
              loggamma(i, k) = calculate_gamma(eta_star_curr.slice(tt),
                       z, beta, curr_clu(tt, i)-1, k, i, 1);
            }
          }else{
            /*id_empty: it identifies which of the augmented has been chosen
             * 2b  page 345 Favaro and Teh, second column
             */
            id_empty = newci - nclu_curr(tt) - 1;
            nclu_curr(tt) += 1;
            curr_clu(tt, i) = nclu_curr(tt);
            nj_curr(tt, nclu_curr(tt)-1) = 1;

            eta_star_curr.slice(tt).col(nclu_curr(tt)-1) = eta_empty.slice(tt).col(id_empty);

            //NEED TO UPDATE GAMMA TOO
            //for(ii = 0; ii < nobs; ii++){
            for(k = 0; k < dim; k++){
              loggamma(i, k) = calculate_gamma(eta_star_curr.slice(tt),
                       z, beta, curr_clu(i)-1, k, i, 1);
            }
            //}
            eta_empty.slice(tt).col(id_empty) = ran_mvnorm(hP0_m0, hP0_L0, dim);
          }

          if(PPMx == 1){
            // need to now add the xcon to the cluster to which it was assigned;
            for(p = 0; p < ncon; p++){
              sumx(tt, (curr_clu(tt, i)-1)*ncon + p) += xcon(i * ncon + p);
              sumx2(tt, (curr_clu(tt, i)-1)*ncon + p) += xcon(i * ncon + p) * xcon(i * ncon + p);
            }
            for(p = 0; p < ncat; p++){
              njc(tt, ((curr_clu(tt, i)-1)*ncat + p) * max_C + xcat(i * ncat + p)) += 1;
            }
          }
        }
      }

      //////////////////////////////////////////////////
      // update the cluster value with NEAL 8 (II step)
      //////////////////////////////////////////////////
      /*Sthetam = (eta_star_curr.col(j) - thetaj)*(eta_star_curr.col(j) - thetaj).t();
       Sthetam = arma::inv(S0m + Sthetam);
       Sigma = ran_iwish(nuiw + 1, Sthetam, dim);*/

      for(j = 0; j < nclu_curr(tt); j++){
        //Rcpp::Rcout << "lg_eu" << loggamma << std::endl;
        l_eta_out = eta_update(JJ, loggamma, nclu_curr(tt), curr_clu.row(tt).t(), nj_curr.row(tt).t(),
                               treatments, tt, eta_star_curr.slice(tt).col(j), eta_flag.row(tt).t(),
                               mu0, L0v, j, mhtunepar(0));
        eta_star_curr.slice(tt).col(j) = Rcpp::as<arma::vec>(l_eta_out[0]);
        loggamma = Rcpp::as<arma::mat>(l_eta_out[1]);
        eta_flag.row(tt) = Rcpp::as<arma::vec>(l_eta_out[2]);
      }
      sumtotclu(tt) += nclu_curr(tt);

      if(upd_hier == 1){
        for(j = 0; j < (dim*dim); j++){
          Swork(j) = 0.0;
        }

        for(j = 0; j < nclu_curr(tt); j++){
          Swork += (eta_star_curr.slice(tt).col(j) - emme0) * (eta_star_curr.slice(tt).col(j) - emme0).t();
        }

        Bwork = nuiw*Psi0_mat + Swork;
        Bwork = arma::inv(Bwork);
        for(i = 0; i < dim; i++){
          for(j = 0; j < dim; j++){
            Bvecwork(i*dim + j) = Bwork(i, j);
          }
        }

        L0v = ran_iwish(nclu_curr(tt) + nuiw, Bvecwork, dim);

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

        for(j = 0; j < nclu_curr(tt); j++){
          Rwork += eta_star_curr.slice(tt).col(j);
        }

        //Rcpp::Rcout << "Rwork" << Rwork.t() << std::endl;

        Awork = Vwork*((arma::inv(sigma0_mat)*emme0)+ arma::inv(L0_mat)*Rwork);

        for(i = 0; i < dim; i++){
          for(j = 0; j < dim; j++){
            Vvecwork(i*dim + j) = Vwork(i, j);
          }
        }

        mu0 = ran_mvnorm(Awork, Vvecwork, dim);
        //Rcpp::Rcout << "mu0" << mu0.t() << std::endl;
      }

      for(k = 0; k < dim; k++){
        for(q = 0; q < Q; q++){
          h = q + k * Q;
          //h = k + q * dim;
          beta_temp(q) = beta(h);
        }

        l_beta_out = beta_update(z, JJ, loggamma, beta_temp, beta_flag,
                                 treatments, tt, mu_beta,
                                 sigma_beta, k, mhtunepar(1));

        beta_temp = Rcpp::as<arma::vec>(l_beta_out[0]);
        for(q = 0; q < Q; q++){
          h = q + k * Q;
          //h = k + q * dim;
          beta(h) = beta_temp(q);
        }
        loggamma = Rcpp::as<arma::mat>(l_beta_out[1]);
        beta_flag = Rcpp::as<arma::mat>(l_beta_out[2]);
      }

      if(hsp == 1){
        lambda = up_lambda_hs(beta, lambda, tau);
        tau = up_tau_hs(beta, lambda, tau);
        sigma_beta = lambda * tau;
      }

      /*////////////////////////////////////////////////
       * update random variables:
       * - independent gammas for DM sampling scheme JJ, TT
       * - random gamma (trucchetto Raffaele) ss
       ////////////////////////////////////////////////*/
       // update JJ and consequently TT
       for(i = 0; i < nobs; i++){
         if(treatments(i) == tt){
           TT(i) = 0.0;
           for(k = 0; k < dim; k++){
             JJ(i, k) = R::rgamma(y(i, k) + exp(loggamma(i, k)), pow(ss(i) + 1.0, - 1.0));
             if(JJ(i, k) < pow(10.0, -100.0)){
               JJ(i, k) = pow(10.0, -100.0);
             }
             TT(i) += JJ(i, k);
           }
         }
       }

       // update latent variables ss
       for(i = 0 ; i < nobs ; i++){
         if(treatments(i) == tt){
           ss(i) = R::rgamma(ypiu(i), pow(TT(i), - 1.0));
         }
       }

       /*////////////////////////////////////////////////
        * in sample prediction to assess model fit
        ////////////////////////////////////////////////*/
        if((l > (burn-1)) & ((l + 1) % thin == 0)){
          if(treatments(i) == tt){
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
        }
    }//this closes the loop on candidate therapies T
     // for the training set

    /*////////////////////////////////////////////////
     * Posterior Predictive
     ////////////////////////////////////////////////*/
     if((l > (burn-1)) & ((l + 1) % thin == 0)){
       for(pp = 0; pp < npred; pp++){//loop for every subject in the test set
         for(tt = 0; tt < nT; tt++){
           for(j = 0; j < nclu_curr(tt); j++){

             lgconN=0.0, lgconY=0.0;
             lgcatN=0.0, lgcatY=0.0;

             if(PPMx == 1){
               for(p = 0; p < ncon; p++){
                 sumxtmp = sumx(tt, j * ncon + p);
                 sumx2tmp = sumx2(tt, j * ncon + p);
                 if(similarity==1){ // Auxilliary
                   if(consim==1){//normal normal
                     lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(tt, p), nj_curr(tt, j), 0, 0, 1);
                     lgconN += lgcont;
                   }
                   if(consim==2){//normal normal inverse gamma
                     lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(tt, p), s2mle(tt, p), nj_curr(tt, j), 0, 0, 1);
                     //sufficient statistics da controllare: in Page "sballate"
                     lgconN += lgcont;
                   }
                 }
                 if(similarity==2){ //Double Dipper
                   if(consim==1){//normal normal
                     lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(tt, p), nj_curr(tt, j), 1, 0, 1);
                     //sufficient statistics da controllare: in Page "sballate"
                     lgconN += lgcont;
                   }
                   if(consim==2){//normal normal inverse gamma
                     lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(tt, p), s2mle(tt, p), nj_curr(tt, j), 1, 0, 1);
                     lgconN += lgcont;
                   }
                 }

                 //add the pp-th predictor to cluster
                 sumxtmp += xconp(pp * ncon + p);
                 sumx2tmp += xconp(pp * ncon + p) * xconp(pp * ncon + p);

                 if(similarity==1){ // Auxilliary
                   if(consim==1){//normal normal
                     lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(tt, p), nj_curr(tt, j) + 1, 0, 0, 1);
                     lgconY += lgcont;
                   }
                   if(consim==2){//normal normal inverse gamma
                     lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(tt, p), s2mle(tt, p), nj_curr(tt, j) + 1, 0, 0, 1);
                     //sufficient statistics da controllare: in Page "sballate"
                     lgconY += lgcont;
                   }
                 }
                 if(similarity==2){ //Double Dipper
                   if(consim==1){//normal normal
                     lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, xbar(tt, p), nj_curr(tt, j) + 1, 1, 0, 1);
                     //sufficient statistics da controllare: in Page "sballate"
                     lgconY += lgcont;
                   }
                   if(consim==2){//normal normal inverse gamma
                     lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, xbar(tt, p), s2mle(tt, p), nj_curr(tt, j) + 1, 1, 0, 1);
                     lgconN += lgcont;
                   }
                 }
               }//this closes the loop on continuous covariates

               for(p = 0; p < ncat; p++){

                 for(c = 0; c < max_C; c++){
                   njctmp(c) = njc(tt, (j*ncat + p)*(max_C) + c);
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

             weight(tt, j) = log((double) nj_curr(tt, j)) + // cohesion part
               lgcatY - lgcatN + // Categorical part
               lgconY - lgconN;  // Continuous part

             if(calibration == 2){
               weight(tt, j) = log((double) nj_curr(tt, j)) + // cohesion part
                 (1/((double)ncon + (double)ncat))*(lgcatY + lgconY - lgcatN - lgconN);
             }
           }//this closes the loop on existing clusters

           lgcondraw = 0.0;
           lgcatdraw = 0.0;

           //probabilità che la predittiva metta l'osservazione in un cluster tutto suo
           for(j = nclu_curr(tt); j < (nclu_curr(tt) + CC); j++){
             jj = j - nclu_curr(tt);
             if(PPMx == 1){
               // Continuous Covariates
               for(p = 0; p < (ncon); p++){
                 tmp = xconp(pp*(ncon) + p);
                 if(similarity==1){ // Auxilliary
                   if(consim==1){//normal normal
                     lgcondraw += gsimconNN(m0, v, s20, tmp, tmp*tmp, xbar(tt, p), 1, 0, 0, 1);
                   }
                   if(consim==2){//normal normal inverse gamma
                     lgcondraw += gsimconNNIG(m0, k0, nu0, s20, tmp, tmp*tmp, xbar(tt, p), s2mle(tt, p), 1, 0, 0, 1);
                   }
                 }
                 if(similarity==2){ //Double Dipper
                   if(consim==1){//normal normal
                     lgcondraw += gsimconNN(m0, v, s20, tmp, tmp*tmp, xbar(tt, p), 1, 1, 0, 1);
                   }
                   if(consim==2){//normal normal inverse gamma
                     lgcont += gsimconNNIG(m0, k0, nu0, s20, tmp, tmp*tmp, xbar(tt, p), s2mle(tt, p), 1, 1, 0, 1);
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

             weight(tt, j) = log(alpha) - log(CC) + //cohesion + auxiliary ptms
               lgcondraw + // Continuous covariate part
               lgcatdraw; // categorical covariate part
             if(calibration == 2){
               weight(tt, j) = log(alpha) - log(CC) +
                 (1/((double)ncon + (double)ncat))*(lgcondraw + lgcatdraw);
             }
           }//chiude loop su empty cluster

           if((calibration == 1) & (PPMx == 1)){
             maxgtilN = gtilN(0);
             maxgtilY = gtilY(0);

             for(j = 1; j < nclu_curr(tt) + CC; j++){

               if(maxgtilN < gtilN(j)) maxgtilN = gtilN(j);

               if(j < nclu_curr(tt)){
                 if(maxgtilY < gtilY(j)) maxgtilY = gtilY(j);
               }
             }

             sgY = 0.0;
             sgN = 0.0;

             for(j = 0; j < nclu_curr(tt) + CC; j++){

               lgtilN(j) = gtilN(j) - maxgtilN;
               sgN = sgN + exp(lgtilN(j));

               if(j < nclu_curr(tt)){// If x is included in an existing cluster in cannot be a singleton
                 lgtilY(j) = gtilY(j) - maxgtilY;
                 sgY = sgY + exp(lgtilY(j));
               }
             }
             // Calibrate the unnormalized cluster probabilities
             for(j = 0; j < nclu_curr(tt); j++){
               lgtilNk = lgtilN(j) - log(sgN);
               lgtilYk = lgtilY(j) - log(sgY);

               weight(tt, j) = log((double) nj_curr(tt, j)) +  // Cohesion part
                 lgtilYk - lgtilNk; //This takes into account both cont and cat vars
             }

             // calibration for empty clusters
             for(j = nclu_curr(tt); j < nclu_curr(tt) + CC; j++){
               jj = j - nclu_curr(tt);
               lgtilNk = lgtilN(j) - log(sgN);
               //lgtilYk = lgtilY(j) - log(sgY);

               weight(tt, j) = log(alpha) - log(CC) +  // Cohesion part
                 //lgtilYk -
                 lgtilNk;
               //Rcpp::Rcout << "changed" << std::endl;
               /*weight(j) = log(alpha) - log(CC) +
                lgtilN(j) - // Continuous covariate part
                log(sgN);*/
             }
           }//chiude calibrazione 1

           //AVOID ZERO IN WEIGHTS
           maxwei = weight(tt, 0);
           for(j = 1; j < nclu_curr(tt) + CC; j++){
             if(maxwei < weight(tt, j)) maxwei = weight(tt, j);
           }

           denwei = 0.0;

           for(j = 0; j < nclu_curr(tt) + CC; j++){
             weight(tt, j) = exp(weight(tt, j) - maxwei);
             denwei += weight(tt, j);
           }

           for(j = 0; j < nclu_curr(tt) + CC; j++){
             pweight(tt, j) = weight(tt, j)/denwei;
             //mysws += pweight(j);
           }

           //sample the new cluster for i-th observation
           uu = R::runif(0.0,1.0);
           cweight = 0.0;
           //newci = id_empty;
           for(j = 0; j < nclu_curr(tt) + CC; j++){
             cweight += pweight(tt, j);

             if (uu < cweight){
               newci = j + 1;
               break;
             }
           }

           /*
            * adjust cluster labels and cardinalities
            */

           if((newci) <= (nclu_curr(tt))){
             eta_pred = eta_star_curr.slice(tt).col(newci - 1);
           }else{
             eta_pred = ran_mvnorm(mu0, L0v, dim);
           }

           //NEED TO UPDATE GAMMA TOO
           for(k = 0; k < dim; k++){
             loggamma_pred(k) = calculate_gamma(eta_pred, zpred, beta, 0, k, pp, 1);
           }

           pii_pred.slice(tt).row(pp) = exp(loggamma_pred).t();

           pii_pred.slice(tt).row(pp) = pii_pred.slice(tt).row(pp)/arma::sum(pii_pred.slice(tt).row(pp));

           ppred.slice(tt).row(pp) = rmultinom_rcpp(1, 1, pii_pred.slice(tt).row(pp).t());
           predclass(pp, tt) = newci;
         }//this closesthe loop for each treatment
       }//this closes the loop on npred subjects
     }//this closes the if for the draws after burnin and thinned

     //////////////////////
     // Save MCMC iterates
     //////////////////////
     if((l > (burn-1)) & ((l + 1) % thin == 0)){
       nclus(ll) = nclu_curr;
       for(i = 0; i < nobs; i++){
         Clui(ll*(nobs) + i) = curr_clu(i);
         like(ll*(nobs) + i) = like_iter(i);
         beta_out.col(ll) = beta;
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
    Rcpp::Named("beta") = beta_out,
    Rcpp::Named("eta_acc") = eta_flag/sumtotclu,
    Rcpp::Named("beta_acc") = beta_flag,
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
