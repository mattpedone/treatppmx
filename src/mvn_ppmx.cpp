#include <RcppArmadillo.h>
#include <math.h>
#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include "utils.h"

// [[Rcpp::export]]
Rcpp::List mvn_ppmx(int iter, int burn, int thin, int nobs, int ncon, int ncat, 
                    arma::vec catvec, double alpha, int CC, int cohesion, 
                    int similarity, arma::mat y, arma::vec xcon, 
                    arma::vec xcat, arma::vec similparam, //arma::vec modelpriors, 
                    arma::vec hP0_m0, arma::vec hP0_L0, double hP0_nu0, 
                    arma::vec hP0_V0,
                    arma::vec mhtune, int calibration){
  
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
  int l, ll, i, ii, c, p, j, mm, zi;
  
  int dim = y.n_cols;
  double max_C, nout, sumx, sumx2;
  double ldet0; //for mvn density
  //number of saved iterations
  nout = (iter - burn)/(thin);
  Rcpp::Rcout << "nout =  " << nout << std::endl;
  //numero totale di covariate
  //int ncov = ncon + ncat;
  
  //maximum number of categories for categorical covariates
  max_C = catvec.max();
  
  //sum of response vector
  arma::vec sumy(dim); 
  sumy = sum(y, 0);
  
  //compute mean and variance of each continuous covariate
  arma::vec mnmle(ncon);
  arma::vec s2mle(ncon);
  for(p = 0; p < ncon; p++){
    sumx = 0.0, sumx2 = 0.0;
    for(i = 0; i < nobs; i ++){
      sumx = sumx + xcon(i*(ncon) + p);
      sumx2 = sumx2 + xcon(i*(ncon) + p)*xcon(i*(ncon) + p);
    }
    
    mnmle(p) = sumx/((double) nobs);
    s2mle(p) = sumx2/((double) nobs) - mnmle(p)*mnmle(p);
  }
  
  //////////////////////////
  ////Cluster-related stuff
  //////////////////////////
  
  //upper bound on sig20 
  //that is upper bound of uniform prior on cluster mean variance
  
  /*inizializzo il vec delle label dei cluster, il vettore delle cardinalità
   * di ciascun cluster e l'int che mi dice il numero di cluster correnti
   */
  arma::vec curr_clu(nobs);
  for(i = 0; i < nobs; i++){
    curr_clu(i) = R::rbinom(2, 0.25) + 1;
  }
  
  arma::vec nj_curr(nobs);
  nj_curr.fill(0);
  //cfr adr1
  for(i = 0; i < nobs; i++){
    for(ii = 0; ii < nobs; ii++){
      if(curr_clu(i) == ii+1) nj_curr(ii) = nj_curr(ii) + 1;
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
  
  ////////////////////////////////////////////////////
  //// cluster specific parameters stuff
  ////////////////////////////////////////////////////
  
  arma::mat mu_star_curr(nclu_curr, dim);
  mu_star_curr.fill(0); //tbfilled from input
  
  arma::mat sigma_star_curr(nclu_curr, dim * dim);
  sigma_star_curr.fill(1); //tbfilled from input
  
  ////////////////////////////////////////////////////
  //// reuse algorithm stuff
  ////////////////////////////////////////////////////
  
  arma::mat mu_empty(CC, dim, arma::fill::zeros);
  arma::mat sigma_empty(CC, dim * dim, arma::fill::ones);
  
  /*double sig20_iter = 0.5 * modelpriors(3); 
   
   
   int iaux;
   
   
   arma::vec njc((nobs)*(ncat));
   
   arma::vec xcontmp(nobs);
   
   
   
   
   arma::vec sig2h(nobs); 
   sig2h.fill(0.5*modelpriors(3));
   arma::vec muh(nobs);
   muh.fill(0.0);
   */
  /////////////////////////////////////
  // Storage for posterior predictive
  /////////////////////////////////////
  
  //arma::vec ispred_iter(nobs); 
  //ispred_iter.fill(0.0);
  
  ////////////////////////////////////////
  // Storage needed to update parameters
  ////////////////////////////////////////
  
  // update Si (cluster labels);
  /*int njtmp;
   double auxm, auxs2, tmp;//, npdN, npdY, npd;
   double mudraw, sdraw, maxph, denph, cprobh, uu;
   double lgconN, lgconY, lgcatN, lgcatY, lgcondraw, lgcatdraw;
   double lgcont, lgcatt;
   */
  
  /*arma::vec muaug(CC);
   muaug.fill(0.0);
   arma::vec saug(CC);
   saug.fill(0.0);
   
   arma::vec ph(nobs + CC);
   ph.fill(0.0);
   arma::vec probh(nobs + CC);
   probh.fill(0.0);
   
   arma::vec gtilN(nobs + CC);
   gtilN.fill(0.0);
   arma::vec gtilY(nobs + CC);
   gtilY.fill(0.0);
   arma::vec lgtilN(nobs + CC);
   lgtilN.fill(0.0);
   arma::vec lgtilY(nobs + CC);
   lgtilY.fill(0.0);
   
   double sgY, sgN,  lgtilNk, lgtilYk, maxgtilY, maxgtilN;
   
   */
  // update muh (mustar) 
  
  /*
   // update mu0 (mean mustar)
   double summu;
   
   // update sig20 (variance of mustar)
   double llo, lln, llr, os0, ns0;
   
   // update sig2 (sigmastar)
   double osig, nsig;
   */
  // Stuff for out of sample predictions 
  //double lgcon0, lgcat0, mupred, sig2pred;
  
  // Stuff to compute lpml (log pseudo marginal likelihood), 
  // likelihood, and WAIC widely applicable information criterion (WAIC), 
  // also known as Watanabe–Akaike information criterion
  /*double lpml_iter, elppdWAIC;
   arma::vec CPOinv(nobs);
   CPOinv.fill(0.0);
   arma::vec like_iter(nobs);
   like_iter.fill(0.0);
   arma::vec mnlike(nobs);
   mnlike.fill(0.0);
   arma::vec mnllike(nobs);
   mnllike.fill(0.0);
   */
  //	Hyper-prior parameters 
  /*
   // priors for mu0
   double m = modelpriors(0); 
   double s2 = modelpriors(1);
   
   // prior values for sig2h and sig20;
   double smin=0, smax = modelpriors(2);
   double s0min=0, s0max = modelpriors(3);
   
   // DP weight parameter
   double alphadp = alpha;
   
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
   
   // For the variance similarity function
   //double alpha_vs = similparam(6);
   
   // M-H tuning parameters
   double csigSIG0 = mhtune(0);
   double csigSIG = mhtune(1);
   */
  
  /*for(mm = 0; mm < CC; mm++){
   muaug(mm) = R::rnorm(mu0_iter, sqrt(sig20_iter));
   saug(mm) = R::runif(smin, smax);
  }*/
  
  //storage for return 
  /*arma::vec mu(nout * nobs, arma::fill::ones); 
   arma::vec sig2(nout * nobs, arma::fill::ones); 
   arma::vec Si(nout * nobs, arma::fill::ones); 
   arma::vec like(nout * nobs, arma::fill::ones); 
   arma::vec ispred(nout * nobs, arma::fill::ones); 
   
   arma::vec mu0(nout, arma::fill::ones); 
   arma::vec sig20(nout, arma::fill::ones);
   arma::vec nclus(nout, arma::fill::ones);
   
   double WAIC = 1.0;
   double lpml = 1.0;
   
   ll = 0;
   */
  ////////////////////////////////////////
  //
  // HERE COMES THE MCMC
  //
  ////////////////////////////////////////  
  
  for(l = 0; l < iter; l++){
    
    if((l+1) % 10000 == 0){
      Rcpp::Rcout << "mcmc iter =  " << l+1 << std::endl;
    }
    
    /* inizializzo latent variables per empty components da P0 
     * per il Reuse algorithm (vettori di lunghezza CC)
     */
    
    for(mm = 0; mm < CC; mm++){
      mu_empty.row(mm) = ran_mvnorm(hP0_m0, hP0_L0, dim);
      sigma_empty.row(mm) = ran_wish(hP0_nu0, hP0_V0, dim);
    }
    
    /////////////////////////////////////////
    // update the cluster labels with NEAL 8 
    /////////////////////////////////////////
    
    for(i = 0; i < nobs; i++){
      
      zi = curr_clu(i)-1; //sottraggo 1 perché C conta da 0
      
      if(nj_curr(zi) > 1){// Case where the nj corresponding to zi is such that nj>1
        /* l'osservazione NON è un singoletto 
         * Neal lo chiama n_{-1, c}, 
         * nel modello di Raffaele S_{\tilde{j}}^{a\star, -i}
         */
        
        nj_curr(zi) -= 1;
        
        /*
         * - nei pesi che calcolo dopo c'è solo la densità, mancano log(nj), 
         * log(alpha), log(CC)
         * - controlla nel file in cui faccio il confronto con Hoff se 
         * la varianza è specificata bene
         * - in particolare la funzione dinvwish non è corretta (cfr chpt7hoff.R)
         * - questione matrici varianza/precision
         * - forse la cosa migliore è scrivere una funzioncina che faccia inversa e usarla x le densità
         * - verifica se ldet0 va bene
         * - bisogna inserire (capisci se direttamente nel ciclo) le covariate
         * attraverso le funzioni di similarità
         */
        ldet0 = -dim*log(hP0_L0(0));
        
        for(j = 0; j < nclu_curr + CC; j++){
          if(j < nclu_curr){
            //vogliono le inverse!!!
            weight(j) = dmvnorm(y.row(i), mu_star_curr.row(j), 
                   sigma_star_curr.row(j), dim, ldet0, 1);
          } else {
            weight(j) = dmvnorm(y.row(i), mu_empty.row(j), 
                   sigma_empty.row(j), dim, ldet0, 1);
          }
        }
        
        /*
         * qui ci va la normalizzazione dei pesi con l'exp per evitare gli 0
         * cfr ll. 68-70 raffaele
         * e il punto 2a pag 345 Favaro and Teh
         */
        
      }else{
        // l'osservazione è un singoletto 
        iaux = curr_clu(i);
        if(iaux < nclu_curr){
          //questa condizione è vera se il singoletto "sta in mezzo" a non-singleton clusters
          //o meglio il singoletto non sta nell'ultimo cluster
          //quindi devo fare il relabel dei clusters
          //
          // faccio lo swap cluster labels tra curr_clu(i) e nclu_curr
          // ANCHE cluster specific parameters;
          
          // Metto tutti i membri dell'ultimo cluster nel cluster del soggetto i (iaux)
          for(ii = 0; ii < nobs; ii++){
            if(curr_clu(ii) == nclu_curr){
              curr_clu(ii) = iaux;
            }
          }
          
          //metto il soggetto i (singoletto) come ultimo cluster
          curr_clu(i) = nclu_curr;
          
          // The following steps swaps order of cluster specific parameters
          // so that the newly labeled subjects from previous step retain
          // their correct cluster specific parameters
          
          auxm = muh(iaux-1);
          muh(iaux-1) = muh(nclu_curr-1);
          muh(nclu_curr-1) = auxm;
          
          auxs2 = sig2h(iaux-1);
          sig2h(iaux-1) = sig2h(nclu_curr-1);
          sig2h(nclu_curr-1) = auxs2;
          
          
          // the number of members in cluster is also swapped with the last
          nj_curr(iaux-1) = nj_curr(nclu_curr-1);
          nj_curr(nclu_curr-1) = 1;
        }
        
        
        // Now remove the ith obs and last cluster;
        nj_curr(nclu_curr-1) = nj_curr(nclu_curr-1) - 1;
        nclu_curr = nclu_curr - 1;
      }
      
      // The atoms have been relabeled 
      // Begin the cluster probabilities
      
      for(j = 0; j < nclu_curr; j++){
        
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
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcont = similarityf::gsimconNN(m0, v, s20, sumx, sumx2, mnmle(p), njtmp, 0, 0, 1);
              lgconN = lgconN + lgcont;
            }
            if(consim==2){//normal normal inverse gamma
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcont = similarityf::gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle(p), s2mle(p), njtmp, 0, 0, 1);
              lgconN = lgconN + lgcont;
            }
            
          }
          if(similarity==2){ //Double Dipper
            if(consim==1){//normal normal
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcont = similarityf::gsimconNN(m0, v, s20, sumx, sumx2, mnmle(p), njtmp, 1, 0, 1);
              lgconN = lgconN + lgcont;
            }
            if(consim==2){//normal normal inverse gamma
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcont = similarityf::gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle(p), s2mle(p), njtmp, 1, 0, 1);
              lgconN = lgconN + lgcont;
            }
          }
          
          // now add ith individual back;
          
          xcontmp(njtmp) = xcon(i*(ncon)+p);
          sumx = sumx + xcon(i*(ncon)+p);
          sumx2 = sumx2 + xcon(i*(ncon)+p)*xcon(i*(ncon)+p);
          njtmp += 1;
          
          if(similarity==1){ // Auxilliary
            if(consim==1){//normal normal
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcont = similarityf::gsimconNN(m0, v, s20, sumx, sumx2, mnmle(p), njtmp, 0, 0, 1);
              lgconY = lgconY + lgcont;
            }
            if(consim==2){//normal normal inverse gamma
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcont = similarityf::gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle(p), s2mle(p), njtmp, 0, 0, 1);
              lgconY = lgconY + lgcont;
            }
          }
          if(similarity==2){ //Double Dipper
            if(consim==1){//normal normal
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcont = similarityf::gsimconNN(m0, v, s20, sumx, sumx2, mnmle(p), njtmp, 1, 0, 1);
              lgconY = lgconY + lgcont;
            }
            if(consim==2){
              //normal normal inverse gamma
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcont = similarityf::gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle(p), s2mle(p), njtmp, 1, 0, 1);
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
          
          
          if(similarity == 1){ // Auxiliary
            //Dirichlet-Multinomial
            //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
            lgcatt = similarityf::gsimcatDM(njc, dirweights, catvec(p), 0, 1);
            lgcatN = lgcatN + lgcatt;
          }
          if(similarity==2){// Double dipper
            //Dirichlet-Multinomial
            //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
            lgcatt = similarityf::gsimcatDM(njc, dirweights, catvec(p), 1, 1);
            lgcatN = lgcatN + lgcatt;
          }
          
          njc(xcat(i*(ncat)+p)) += 1;
          njtmp += 1;
          
          if(similarity==1){
            //Dirichlet-Multinomial
            //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
            lgcatt = similarityf::gsimcatDM(njc, dirweights, catvec(p), 0, 1);
            lgcatY = lgcatY + lgcatt;
          }
          if(similarity==2){
            //Dirichlet-Multinomial
            //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
            lgcatt = similarityf::gsimcatDM(njc, dirweights, catvec(p), 1, 1);
            lgcatY = lgcatY + lgcatt;
          }
        }//chiude ciclo su p covariate discrete
        
        gtilY(j) = lgconY + lgcatY;
        gtilN(j) = lgconN + lgcatN;
        
        //ASTERISCO
        ph(j) = R::dnorm(y(i), muh(j), sqrt(sig2h(j)), 1) +
          log((double) nj_curr(j)) + // cohesion part
          lgcatY - lgcatN + // Categorical part
          lgconY - lgconN;  // Continuous part
        
        if(calibration == 2){
          ph(j) = R::dnorm(y(i), muh(j), sqrt(sig2h(j)), 1) +
            log((double) nj_curr(j)) + // cohesion part
            (1/((double)ncon + (double)ncat))*(lgcatY + lgconY - lgcatN - lgconN);
        }
        
        // Uniform cohesion
        if(cohesion==2){
          ph(j) = ph(j) - log((double) nj_curr(j));
        }
      }//chiude il ciclo sui j cluster
      
      // Now need to generate auxiliary variable from the prior so that there is
      // possibility of new cluster
      
      //EXT qui campioni m valori
      //mudraw e sdraw li usi alla fine per salvare il valore eventualmente effettivamente
      //pescato dagli augmented auxiliary parameters
      
      for(mm = 0; mm < CC; mm++){
        muaug(mm) = R::rnorm(mu0_iter, sqrt(sig20_iter));
        saug(mm) = R::runif(smin, smax);
      }
      
      //qui non dovrebbe esserci niente da modificare
      // Continuous Covariates
      lgcondraw = 0.0;
      for(p = 0; p < (ncon); p++){
        xcontmp(0) = xcon(i*(ncon)+p);
        
        if(similarity == 1){ // Auxiliary
          if(consim == 1){
            lgcont = similarityf::gsimconNN(m0,v,s20,xcontmp(0),xcontmp(0)*xcontmp(0), mnmle(p),1,0,0, 1);
            lgcondraw = lgcondraw + lgcont;
          }
          if(consim == 2){
            lgcont = similarityf::gsimconNNIG(m0, k0, nu0, s20, xcontmp(0), xcontmp(0)*xcontmp(0),mnmle(p),s2mle(p), 1, 0,0, 1);
            lgcondraw = lgcondraw + lgcont;
          }
        }
        if(similarity == 2){ // Double Dipper
          if(consim == 1){
            lgcont = similarityf::gsimconNN(m0,v,s20,xcontmp(0),xcontmp(0)*xcontmp(0), mnmle(p), 1, 1, 0, 1);
            lgcondraw = lgcondraw + lgcont;
          }
          if(consim == 2){
            lgcont = similarityf::gsimconNNIG(m0, k0, nu0, s20, xcontmp(0), xcontmp(0)*xcontmp(0),mnmle(p),s2mle(p), 1, 1, 0, 1);
            lgcondraw = lgcondraw + lgcont;
          }
        }
      }
      
      //qui non dovrebbe esserci niente da modificare
      // Categorical Covariates
      lgcatdraw = 0.0;
      for(p = 0; p < (ncat); p++){
        for(c = 0; c < catvec(p); c++){
          njc(c) = 0;
        }
        
        njc(xcat(i*(ncat)+p)) = 1;
        
        if(similarity == 1){
          lgcatt = similarityf::gsimcatDM(njc, dirweights, catvec(p), 0, 1);
          lgcatdraw = lgcatdraw + lgcatt;
        }
        if(similarity == 2){
          lgcatt = similarityf::gsimcatDM(njc, dirweights, catvec(p), 1, 1);
          lgcatdraw = lgcatdraw + lgcatt;
        }
      }
      
      //qui metti ciclo dopo cicli su p ncon e p ncat
      for(mm = nclu_curr; mm < (nclu_curr+CC); mm++){
        gtilY(mm) = lgcondraw + lgcatdraw;
        gtilN(mm) = lgcondraw + lgcatdraw;
      }
      //gtilY(nclu_curr) = lgcondraw + lgcatdraw;
      //gtilN(nclu_curr) = lgcondraw + lgcatdraw;
      
      //EXT
      for(mm = nclu_curr; mm < (nclu_curr+CC); mm++){
        ph(mm) = R::dnorm(y(i), muaug(mm - nclu_curr), saug(mm - nclu_curr), 1) +
          log(alphadp) +  // DP part
          lgcondraw + // Continuous covariate part
          lgcatdraw; // categorical covariate part
        if(calibration == 2){
          ph(mm) = R::dnorm(y(i), muaug(mm - nclu_curr), saug(mm - nclu_curr), 1) +
            log(alphadp) +
            (1/((double)ncon + (double)ncat))*(lgcondraw + lgcatdraw);
        }
        
        if(cohesion == 2){
          ph(mm) -= log(alphadp);
        }
      }
      
      if(calibration == 1){
        maxgtilN = gtilN(0);//arma::max(gtilN)
        maxgtilY = gtilY(0);
        for(j = 1; j < nclu_curr + CC; j++){//1 - CC
          
          if(maxgtilN < gtilN(j)) maxgtilN = gtilN(j);
          
          if(j < nclu_curr){
            if(maxgtilY < gtilY(j)) maxgtilY = gtilY(j);
          }
        }
        
        sgY=0.0;
        sgN=0.0;
        
        for(j = 0; j < nclu_curr + CC; j++){//1 - CC
          
          lgtilN(j) = gtilN(j) - maxgtilN;
          sgN += exp(lgtilN(j));
          
          if(j < nclu_curr){// If x is included in an existing cluster in cannot be a singleton
            lgtilY(j) = gtilY(j) - maxgtilY;
            sgY = sgY + exp(lgtilY(j));
          }
        }
        
        //qui non modificare MA CONTROLLA!!
        // Calibrate the unnormalized cluster probabilities
        
        for(j = 0; j < nclu_curr; j++){
          lgtilNk = lgtilN(j) - log(sgN);
          lgtilYk = lgtilY(j) - log(sgY);
          
          ph(j) = R::dnorm(y(i), muh(j), sqrt(sig2h(j)), 1) +
            log((double) nj_curr(j)) +  // Cohesion part
            lgtilYk - lgtilNk; //This takes into account both cont and cat vars
          
          if(cohesion == 2){
            ph(j) = ph(j) - log((double) nj_curr(j));
          }
        }
        
        // calibration for a singleton
        // EXT
        for(mm = nclu_curr; mm < (nclu_curr+CC); mm++){
          ph(mm) = R::dnorm(y(i), muaug(mm - nclu_curr), saug(mm - nclu_curr), 1) +
            log(alphadp) +
            lgtilN(mm) - log(sgN) - log(CC);//in caso togli -log(CC)
          
          if(cohesion==2){// Note with a uniform cohesion, for a new cluster
            // the value of log(c({nclu_curr}}) = log(1) = 0;
            ph(mm) -= log(alphadp);
          }
        }
        
      }
      
      //NORMALIZZAZIONE PROBABILITà
      maxph = ph(0);
      for(j = 1; j < nclu_curr + CC; j++){//1 - CC
        if(maxph < ph(j)) maxph = ph(j);
      }
      
      denph = 0.0;
      for(j = 0; j < nclu_curr + CC; j++){//1 - CC
        ph(j) = exp(ph(j) - maxph);
        denph += ph(j);
      }
      
      for(j = 0; j < nclu_curr + CC; j++){//1 - CC
        probh(j) = ph(j)/denph;
      }
      
      uu = R::runif(0.0,1.0);
      
      //visto che al massimo posso aggiungere un solo cluster alla volta 
      //lascerei +1 invece che + CC
      cprobh= 0.0;
      iaux = nclu_curr + 1;
      for(j = 0; j < nclu_curr + 1; j++){
        cprobh = cprobh + probh(j);
        if (uu < cprobh){
          iaux = j + 1;//arbitrariamente prendo il primo cfr Neal(2000)
          break;
        }
      }
      
      if(iaux <= nclu_curr){
        
        curr_clu(i) = iaux;
        nj_curr(curr_clu(i)-1) += 1;
        
        //Rcpp::Rcout << "muaug(0) " << muaug(0) << std::endl;
        //Rcpp::Rcout << "saug(0) " << saug(0) << std::endl;
      } else {
        //qui dal vettore di m auxiliary variables salva mudraw & sdraw
        mudraw = muaug(0);
        sdraw = saug(0);
        //Rcpp::Rcout << "here! " << std::endl;
        
        //Rcpp::Rcout << "here 3! " << std::endl;
        nclu_curr += 1;
        curr_clu(i) = nclu_curr;
        nj_curr(curr_clu(i)-1) = 1;
        
        muh(curr_clu(i)-1) = mudraw;
        sig2h(curr_clu(i)-1) = sdraw*sdraw;
        
      }
      
      // Compute the CPO and lpml using the mixture
      
      like_iter(i) = R::dnorm(y(i), muh[curr_clu(i)-1], sig2h(curr_clu(i)-1), 0);
      
      if((l > (burn-1)) & (l % (thin) == 0)){
        
        // These are needed for WAIC
        mnlike(i) = mnlike(i) + (like_iter(i))/(double) nout;
        mnllike(i) = mnllike(i) + log(like_iter(i))/(double) nout;
        
        CPOinv(i) = CPOinv(i) + (1/(double) nout)*(1/like_iter(i));
      }
    }
    
    //////////////////////////////////////////////////
    // update the cluster value with NEAL 8 (II step)
    //////////////////////////////////////////////////
    
    //muh cluster specific mean
    
    for(j = 0; j < nclu_curr; j++){
      sumy = 0.0;
      for(i = 0; i < nobs; i++){
        if(curr_clu(i) == j+1){
          sumy += y(i);
        }
      }
      
      s2star = 1/((double) nj_curr(j)/sig2h(j) + 1/sig20_iter);
      mstar = s2star*( (1/sig2h(j))*sumy + (1/sig20_iter)*mu0_iter);
      
      
      muh(j) = R::rnorm(mstar, sqrt(s2star));
    }
    
    //mu0 prior mean of muh
    
    summu = 0.0;
    for(j = 0; j < nclu_curr; j++){
      summu = summu + muh(j);
    }
    
    s2star = 1/(((double) nclu_curr/sig20_iter) + (1/s2));
    mstar = s2star*((1/sig20_iter)*summu + (1/s2)*m);
    
    mu0_iter = R::rnorm(mstar, sqrt(s2star));
    
    //sigma0 prior variance for muh
    os0 = sqrt(sig20_iter);
    ns0 = R::rnorm(os0,csigSIG0);
    
    if(ns0 > 0){
      
      lln = 0.0;
      llo = 0.0;
      for(j = 0; j < nclu_curr; j++){
        
        llo = llo + R::dnorm(muh(j), mu0_iter, os0,1);
        lln = lln + R::dnorm(muh(j), mu0_iter, ns0,1);
      }
      
      llo = llo + R::dunif(os0, s0min, s0max, 1);
      lln = lln + R::dunif(ns0, s0min, s0max, 1);
      
      
      llr = lln - llo;
      uu = R::runif(0,1);
      
      if(log(uu) < llr){
        sig20_iter = ns0*ns0;
      }
    }
    
    // Update sig2h  cluster specific variance parameters with metropolis and Uniform prior.
    
    
    for(j = 0; j < nclu_curr; j++){
      osig = sqrt(sig2h(j));
      nsig = R::rnorm(osig,csigSIG);
      if((nsig > 0) & (nsig < smax)){
        
        lln = 0.0;
        llo = 0.0;
        for(i = 0; i < nobs; i++){
          if(curr_clu(i) == j+1){
            llo = llo + R::dnorm(y(i), muh(j), osig,1);
            lln = lln + R::dnorm(y(i), muh(j), nsig,1);
          }
        }
        
        llo = llo + R::dunif(osig, smin, smax, 1);
        lln = lln + R::dunif(nsig, smin, smax, 1);
        
        llr = lln - llo;
        uu = R::runif(0,1);
        
        if(log(uu) < llr){
          sig2h(j) = nsig*nsig;
        }
      }
      
    }
    
    // in sample prediction to assess model fit
    
    if((l > (burn-1)) & (l % (thin) == 0)){
      for(i = 0; i < nobs; i++){
        ispred_iter(i) = R::rnorm(muh(curr_clu(i)-1), sqrt(sig2h(curr_clu(i)-1)));
      }
    }
    
    //////////////////////
    // Save MCMC iterates	
    //////////////////////
    
    if((l > (burn-1)) & ((l + 1) % thin == 0)){
      mu0(ll) = mu0_iter;
      sig20(ll) = sig20_iter;
      
      nclus(ll) = nclu_curr;
      
      for(i = 0; i < nobs; i++){
        mu(ll*(nobs) + i) = muh(curr_clu(i)-1);
        sig2(ll*(nobs) + i) = sig2h(curr_clu(i)-1);
        Si(ll*(nobs) + i) = curr_clu(i);
        
        like(ll*(nobs) + i) = like_iter(i);
        ispred(ll*(nobs) + i) = ispred_iter(i);
      }
      
      ll += 1;
    }  
    
  }//CLOSES MCMC iterations
  
  // calculate LPML
  lpml_iter=0.0;
  
  for(i = 0; i < nobs; i++){
    lpml_iter += log(1/CPOinv(i));
  }
  lpml = lpml_iter;
  
  
  // Computing WAIC  (see Gelman article)
  
  elppdWAIC = 0.0;
  for(i = 0; i < nobs; i++){
    elppdWAIC += (2*mnllike(i) - log(mnlike(i)));
  }
  
  WAIC = -2*elppdWAIC;
  
  //RETURN
  return Rcpp::List::create(Rcpp::Named("mu") = mu,
                            Rcpp::Named("sig2") = sig2,
                            Rcpp::Named("Si") = Si,
                            Rcpp::Named("like") = like,
                            Rcpp::Named("ispred") = ispred,
                            Rcpp::Named("mu0") = mu0,
                            Rcpp::Named("sig20") = sig20,
                            Rcpp::Named("nclus") = nclus,
                            Rcpp::Named("WAIC") = WAIC,
                            Rcpp::Named("lpml") = lpml);
}