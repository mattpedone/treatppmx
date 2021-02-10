#include <RcppArmadillo.h>
#include <math.h>
#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include "utils.h"

//
// [[Rcpp::export]]
Rcpp::List mvn_ppmx(int iter, int burn, int thin, int nobs, int ncon, int ncat,
                    arma::vec catvec, double alpha, int CC, int consim, 
                    int similarity, int calibration, arma::mat y, arma::vec xcon,
                    arma::vec xcat, arma::vec similparam, //arma::vec modelpriors,
                    arma::vec hP0_m0, arma::vec hP0_L0, double hP0_nu0,
                    arma::vec hP0_V0,
                    arma::vec mhtune){
  
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
  double max_C, nout, sumx, sumx2, dval;
  //number of saved iterations
  nout = (iter - burn)/(thin);
  Rcpp::Rcout << "nout =  " << nout << std::endl;
  //numero totale di covariate
  //int ncov = ncon + ncat;

  //maximum number of categories for categorical covariates
  max_C = catvec.max();

  //sum of response vector
  arma::vec sumy(dim);
  sumy = arma::sum(y, 0).t();

  //compute mean and variance of each continuous covariate
  arma::vec xbar(ncon);
  arma::vec s2mle(ncon);
  //arma::vec sumx(ncon);//total of covariates 
  /*arma::vec ssq(ncon * ncon); //matrix stored in an array whose value is \sum_{i=1}^{n} (x_i - xbar) (x_i - xbar)'
   arma::vec sumx2(ncon);*/ 
  
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

  //upper bound on sig20
  //that is upper bound of uniform prior on cluster mean variance

  /*inizializzo il vec delle label dei cluster, il vettore delle cardinalità
   * di ciascun cluster e l'int che mi dice il numero di cluster correnti
   */
  arma::vec curr_clu(nobs);
  for(i = 0; i < nobs; i++){
    curr_clu(i) = R::rbinom(2, 0.25) + 1;
  }
  Rcpp::Rcout << "curr_clu: " << curr_clu << std::endl;

  arma::vec nj_curr(nobs);
  nj_curr.fill(0);
  //cfr adr1
  for(i = 0; i < nobs; i++){
    for(ii = 0; ii < nobs; ii++){
      if(curr_clu(i) == ii+1) nj_curr(ii) += 1;
    }
  }
  Rcpp::Rcout << "nj_curr: " << nj_curr << std::endl;

  int nclu_curr = 0;
  for(i = 0; i < nobs; i++){
    if(nj_curr(i) > 0) nclu_curr += 1;
  }
  Rcpp::Rcout << "nclu_curr: " << nclu_curr << std::endl;
  /*
   * questo è il vettore dei pesi, le probabilità per \rho (Alg 8, Neal 2000)
   * Ha lunghezza pari al massimo numeri di componenti possibile + il numero
   * di variabili ausiliarie che alloco permanentemente con il Reuse Alg
   */

  arma::vec weight(nobs + CC);

  ////////////////////////////////////////////////////
  //// cluster specific parameters stuff
  ////////////////////////////////////////////////////
   
  arma::mat mu_star_curr(nobs, dim);
  mu_star_curr.fill(0); //tbfilled from input

  arma::mat sigma_star_curr(nobs, dim * dim);
  sigma_star_curr.fill(0);
  int idx = 0;
  for(i = 0; i < dim; i++){
    sigma_star_curr(0, idx) = 1;
    idx += (dim + 1);
  }
  for(i = 1; i < nobs; i++){
    sigma_star_curr.row(i) = sigma_star_curr.row(0);
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

  arma::mat mu_empty(CC, dim, arma::fill::zeros);
  arma::mat sigma_empty(CC, dim * dim, arma::fill::zeros);
  idx = 0;
  for(i = 0; i < dim; i++){
    sigma_empty(0, idx) = 1;
    idx += (dim + 1);
  }
  for(i = 1; i < CC; i++){
    sigma_empty.row(i) = sigma_empty.row(0);
  }

  /////////////////////////////////////
  // Storage for posterior predictive
  /////////////////////////////////////
  arma::mat ispred_iter(nobs, dim);
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
  
  // Stuff to compute lpml (log pseudo marginal likelihood),
  // likelihood, and WAIC widely applicable information criterion (WAIC),
  // also known as Watanabe–Akaike information criterion
  double lpml_iter, elppdWAIC;
   arma::vec CPOinv(nobs);
   CPOinv.fill(0.0);
   arma::vec like_iter(nobs);
   like_iter.fill(0.0);
   arma::vec mnlike(nobs);
   mnlike.fill(0.0);
   arma::vec mnllike(nobs);
   mnllike.fill(0.0);

  //Stuff for storage and return
  arma::vec nclus(nout);
  arma::mat mu_out(nout*nobs, dim);
  arma::mat sigma_out(nout*nobs, dim*dim);
  arma::vec Clui(nout * nobs, arma::fill::ones);
  arma::vec like(nout * nobs, arma::fill::ones);
  arma::mat ispred(nout * nobs, dim, arma::fill::ones);
  ll = 0;
  ////////////////////////////////////////
  //
  // HERE COMES THE MCMC
  //
  ////////////////////////////////////////

  for(l = 0; l < iter; l++){
    Rcpp::Rcout << "iterazione" << l+1 << std::endl;
    Rcpp::Rcout << "curr_clu: " << curr_clu << "nj_curr: " << nj_curr << std::endl;
    Rcpp::Rcout << "numero cluster: " << nclu_curr << std::endl;
    if((l+1) % 10000 == 0){
      Rcpp::Rcout << "mcmc iter =  " << l+1 << std::endl;
    }
  
    /* inizializzo latent variables per empty components da P0
     * per il Reuse algorithm (vettori di lunghezza CC)
     */
    
    for(mm = 0; mm < CC; mm++){
      mu_empty.row(mm) = ran_mvnorm(hP0_m0, hP0_L0, dim).t();
      sigma_empty.row(mm) = ran_iwish(hP0_nu0, hP0_V0, dim).t();
    }
    
    for(i = 0; i < nobs; i++){
     // Rcpp::Rcout << "dentro ciclo sulle i" << std::endl;
      /////////////////////////////////////////
      // update the cluster labels with NEAL 8
      /////////////////////////////////////////
      
      zi = curr_clu(i)-1; //sottraggo 1 perché C conta da 0
      //Rcpp::Rcout << "zi" << zi << std::endl;
      if(nj_curr(zi) > 1){// Case where the nj corresponding to zi is such that nj>1
        /* 
         * l'osservazione NON è un singoletto
         * Neal lo chiama n_{-1, c},
         * nel modello di Raffaele S_{\tilde{j}}^{a\star, -i}
         */
        //Rcpp::Rcout << nj_curr(zi) << std::endl;
        nj_curr(zi) -= 1;
        
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
                lgcont = gsimconNN(m0, v, s20, sumx, sumx2, xbar(p), njtmp, 0, 0, 1);
                lgconN += lgcont;
              } 
              if(consim==2){//normal normal inverse gamma
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, xbar(p), s2mle(p), njtmp, 0, 0, 1);
                lgconN += lgcont;
              } 
              
            } 
            if(similarity==2){ //Double Dipper
              if(consim==1){//normal normal
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
                lgcont = gsimconNN(m0, v, s20, sumx, sumx2, xbar(p), njtmp, 1, 0, 1);
                lgconN += lgcont;
              } 
              if(consim==2){//normal normal inverse gamma
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
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
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
                lgcont = gsimconNN(m0, v, s20, sumx, sumx2, xbar(p), njtmp, 0, 0, 1);
                lgconY = lgconY + lgcont;
              } 
              if(consim==2){//normal normal inverse gamma
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, xbar(p), s2mle(p), njtmp, 0, 0, 1);
                lgconY = lgconY + lgcont;
              }
            } 
            if(similarity==2){ //Double Dipper
              if(consim==1){//normal normal
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
                lgcont = gsimconNN(m0, v, s20, sumx, sumx2, xbar(p), njtmp, 1, 0, 1);
                lgconY = lgconY + lgcont;
              } 
              if(consim==2){
                //normal normal inverse gamma
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
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
            
            
            if(similarity == 1){ // Auxiliary
              //Dirichlet-Multinomial
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 0, 1);
              lgcatN = lgcatN + lgcatt;
            } 
            if(similarity==2){// Double dipper
              //Dirichlet-Multinomial
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 1, 1);
              lgcatN = lgcatN + lgcatt;
            } 
            
            njc(xcat(i*(ncat)+p)) += 1;
            njtmp += 1;
            
            if(similarity==1){
              //Dirichlet-Multinomial
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 0, 1);
              lgcatY = lgcatY + lgcatt;
            } 
            if(similarity==2){
              //Dirichlet-Multinomial
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 1, 1);
              lgcatY = lgcatY + lgcatt;
            } 
          }//chiude ciclo su p covariate discrete
          
          gtilY(j) = lgconY + lgcatY;
          gtilN(j) = lgconN + lgcatN;

          /*
           * ogni volta controlla se i determinanti e le matrici di covarianza/precision le passo bene
           */
          
          if(calibration == 1){
            ldSig = logdet(sigma_star_curr.row(j).t(), dim);
            
            weight(j) = /*dmvnorm(y.row(i).t(), mu_star_curr.row(j).t(),
                     sigma_star_curr.row(j).t(), dim, ldSig, 1) +*/
                       log((double) nj_curr(j)) + // cohesion part
                       lgcatY - lgcatN + // Categorical part
                       lgconY - lgconN;  // Continuous part
          }
          if(calibration == 2){
              
            ldSig = logdet(sigma_star_curr.row(j).t(), dim);
            weight(j) = dmvnorm(y.row(i).t(), mu_star_curr.row(j).t(),
                     sigma_star_curr.row(j).t(), dim, ldSig, 1)+
                       log((double) nj_curr(j)) + // cohesion part
                       (1/((double)ncon + (double)ncat))*(lgcatY + lgconY - lgcatN - lgconN);
            }
        }
        
        /*
         * qui ci va la normalizzazione dei pesi con l'exp per evitare gli 0
         * cfr ll. 68-70 raffaele
         * e il punto 2a pag 345 Favaro and Teh
         */
        
        for(j = nclu_curr; j < CC; j++){
          
          // Continuous Covariates
          lgcondraw = 0.0;
          //lgconN = 0.0;
          for(p = 0; p < (ncon); p++){
            tmp = xcon(i*(ncon) + p);
            
          if(similarity==1){ // Auxilliary
              if(consim==1){//normal normal
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
                lgcont = gsimconNN(m0, v, s20, tmp, tmp*tmp, xbar(p), 1, 0, 0, 1);
                lgcondraw += lgcont;
              } 
              if(consim==2){//normal normal inverse gamma
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
                lgcont = gsimconNNIG(m0, k0, nu0, s20, tmp, tmp*tmp, xbar(p), s2mle(p), 1, 0, 0, 1);
                lgcondraw += lgcont;
              } 
              
            } 
          if(similarity==2){ //Double Dipper
              if(consim==1){//normal normal
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
                lgcont = gsimconNN(m0, v, s20, tmp, tmp*tmp, xbar(p), 1, 1, 0, 1);
                lgcondraw += lgcont;
              } 
              if(consim==2){//normal normal inverse gamma
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
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
            
            if(similarity == 1){ // Auxiliary
              //Dirichlet-Multinomial
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 0, 1);
              lgcatdraw += lgcatt;
            } 
            if(similarity==2){// Double dipper
              //Dirichlet-Multinomial
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 1, 1);
              lgcatdraw += lgcatt;
            }
            }//chiude ciclo su covariate discrete
            
          gtilY(j) = lgcondraw + lgcatdraw;
          gtilN(j) = lgcondraw + lgcatdraw;
          /* ogni volta controlla se i determinanti e le matrici di covarianza/precision le passo bene
           */
          
          if(calibration == 1){
            ldSig = logdet(sigma_empty.row(j - nclu_curr).t(), dim);
            weight(j) = dmvnorm(y.row(i).t(), mu_empty.row(j - nclu_curr).t(),
                     sigma_empty.row(j - nclu_curr).t(), dim, ldSig, 1) +
                       log(alpha) - log(CC) +
                       lgcondraw + // Continuous covariate part
                       lgcatdraw; // categorical covariate part
            } 
          if(calibration == 2){
            ldSig = logdet(sigma_empty.row(j - nclu_curr).t(), dim);
            weight(j) = dmvnorm(y.row(i).t(), mu_empty.row(j - nclu_curr).t(),
                     sigma_empty.row(j - nclu_curr).t(), dim, ldSig, 1) +
                       log(alpha) - log(CC) +
                       (1/((double)ncon + (double)ncat))*(lgcondraw + lgcatdraw);
          }
        }
        
        if(calibration == 1){
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
            
            ldSig = logdet(sigma_star_curr.row(j).t(), dim);
            
            weight(j) = dmvnorm(y.row(i).t(), mu_star_curr.row(j).t(),
                   sigma_star_curr.row(j).t(), dim, ldSig, 1) + 
                     log((double) nj_curr(j)) +  // Cohesion part
                     lgtilYk - lgtilNk; //This takes into account both cont and cat vars
            
            /*if(*cohesion==2){
              ph[k] = ph[k] - log((double) nh[k]);
            }*/
          }
          //Rcpp::Rcout << "644" << std::endl;
          // calibration for a singleton
          for(j = nclu_curr; j < nclu_curr + CC; j++){
            /*Rcpp::Rcout << "nclu_curr" << nclu_curr << std::endl;
            int myidx = j- nclu_curr;
            Rcpp::Rcout << "row" << myidx << std::endl;*/
            ldSig = logdet(sigma_empty.row(j - nclu_curr).t(), dim);
            weight(j) = dmvnorm(y.row(i).t(), mu_empty.row(j - nclu_curr).t(),
                   sigma_empty.row(j - nclu_curr).t(), dim, ldSig, 1) +
                     log(alpha) - log(CC) +
                     lgtilN(j) - // Continuous covariate part
                     log(sgN);
            }
          
          /*if(*cohesion==2){// Note with a uniform cohesion, for a new cluster
            // the value of log(c({nclus_iter}}) = log(1) = 0;
            ph[nclus_iter] = ph[nclus_iter] - log(Mdp);
          }*/
          
        }
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
          weight(j) = weight(j)/denwei;
        }
        uu = R::runif(0.0,1.0);
        cweight = 0.0;
        for(j = 0; j < nclu_curr + CC; j++){
          cweight += weight(j);
          
          if (uu < cweight){
            newci = j;
            Rcpp::Rcout << "newci: " << newci << std::endl;
            break;
          }
        }
        if((newci)<nclu_curr){//<=
          //Point 2a page 345 Favaro and Teh, second column
          //Rcpp::Rcout << "RIGA 578ss prima " << std::endl;
          //Rcpp::Rcout << "curr_clu: " << curr_clu << "nj_curr: " << nj_curr << std::endl;
          //Rcpp::Rcout << "numero cluster: " << nclu_curr << std::endl;
          curr_clu(i) = newci + 1;
          nj_curr(newci) += 1;
          //Rcpp::Rcout << "RIGA 579ss dopo " << std::endl;
          //Rcpp::Rcout << "curr_clu: " << curr_clu << "nj_curr: " << nj_curr << std::endl;
          //Rcpp::Rcout << "numero cluster: " << nclu_curr << std::endl;
        }else{
          
          id_empty = newci - nclu_curr;
          //### Point 2b  page 345 Favaro and Teh, second column
          curr_clu(i) = nclu_curr + 1;
          nclu_curr += 1;
          nj_curr(nclu_curr-1) = 1;//nj_curr(nclu_curr-1) = 1;
          //Rcpp::Rcout << "RIGA 589ss " << std::endl;
          //Rcpp::Rcout << "curr_clu: " << curr_clu << "nj_curr: " << nj_curr << std::endl;
          //Rcpp::Rcout << "numero cluster: " << nclu_curr << std::endl;
         /*lascio così perché C conta da 0
          * io voglio mettere 1 nella posizione nclu_curr+1
          * c(nj_current,1);*/
          //## the index of the empty component we have selected is
          //### update the phi parameters
          //Rcpp::Rcout << "punto 598 " << std::endl;
          sigma_star_curr.row(nclu_curr-1) = sigma_empty.row(id_empty);
          //Rcpp::Rcout << "punto 600 " << std::endl;
          mu_star_curr.row(nclu_curr-1) = mu_empty.row(id_empty);
          //Rcpp::Rcout << "punto 602 " << std::endl;
          //## update empty parameter
          mu_empty.row(id_empty) = ran_mvnorm(hP0_m0, hP0_L0, dim).t();
          sigma_empty.row(id_empty) = ran_iwish(hP0_nu0, hP0_V0, dim).t();
        }
      }/* 
        * End of the case we are in a component with nj>1
        */
      else{
        // l'osservazione è un singoletto
        dval = R::runif(0.0, 1.0);
        dval *= (CC);
        id_empty = floor(dval);//praticamente ho fatto un sample()
        
        sigma_empty.row(id_empty) = sigma_star_curr.row(zi);
        mu_empty.row(id_empty) = mu_star_curr.row(zi);
        Rcpp::Rcout << "RIGA 616 " << std::endl;
        Rcpp::Rcout << "curr_clu: " << curr_clu << "nj_curr: " << nj_curr << std::endl;
        Rcpp::Rcout << "numero cluster: " << nclu_curr << std::endl;
        Rcpp::Rcout << "zi: " << zi << std::endl;
        //THIS MAY BE PROBLEMATICO
        for(ii = 0; ii < nobs-1; ii++){//for(ii = 0; ii <= (nobs-1); ii++){
          if(ii >= (zi)){
            nj_curr(ii) = nj_curr(ii+1);
            sigma_star_curr.row(ii) = sigma_star_curr.row(ii + 1);
            mu_star_curr.row(ii) = mu_star_curr.row(ii + 1);
          } //else{nj_curr(ii) = nj_curr(ii);}
        }
        
        nclu_curr -= 1;
        Rcpp::Rcout << "RIGA 635 " << std::endl;
        Rcpp::Rcout << "curr_clu: " << curr_clu << "nj_curr: " << nj_curr << std::endl;
        Rcpp::Rcout << "numero cluster: " << nclu_curr << std::endl;
        Rcpp::Rcout << "zi: " << zi << std::endl;
        
        for(ii = 0; ii < nobs; ii++){
          if(curr_clu(ii)>(zi)){///???
            curr_clu(ii) -= 1;
            /*for(int iii = 0; iii < nclu_curr; iii++){//for(ii = 0; ii <= (nobs-1); ii++){
              if(iii >= (zi)){
                nj_curr(iii) = nj_curr(iii+1);
              }
            }*/
          }
          }
        
        Rcpp::Rcout << "RIGA 650 " << std::endl;
        Rcpp::Rcout << "curr_clu: " << curr_clu << "nj_curr: " << nj_curr << std::endl;
        Rcpp::Rcout << "numero cluster: " << nclu_curr << std::endl;
        
        /*for(ii = 0; ii < nclu_curr; ii++){
          if(ii > zi){
            nj_curr(zi) = nj_curr(zi + 1);
          }
        }*/
        
        /*nj_current <- nj_current[-zi]
        tau_star_current <- tau_star_current[-zi]
        mu_star_current <- mu_star_current[-zi]
        k_current <- k_current-1
        for(ii in 1:n){
          if(clu_current[ii]>zi){
            nj_current <- nj_current[-zi]
          }
        }*/
        //FINE PARTA PROBLEMATICA
        
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
                lgcont = gsimconNN(m0, v, s20, sumx, sumx2, xbar(p), njtmp, 0, 0, 1);
                lgconN += lgcont;
              } 
              if(consim==2){//normal normal inverse gamma
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, xbar(p), s2mle(p), njtmp, 0, 0, 1);
                lgconN += lgcont;
              } 
            
            } 
            if(similarity==2){ //Double Dipper
              if(consim==1){//normal normal
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
                lgcont = gsimconNN(m0, v, s20, sumx, sumx2, xbar(p), njtmp, 1, 0, 1);
                lgconN += lgcont;
              } 
              if(consim==2){//normal normal inverse gamma
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
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
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
                lgcont = gsimconNN(m0, v, s20, sumx, sumx2, xbar(p), njtmp, 0, 0, 1);
                lgconY = lgconY + lgcont;
              } 
              if(consim==2){//normal normal inverse gamma
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, xbar(p), s2mle(p), njtmp, 0, 0, 1);
                lgconY = lgconY + lgcont;
              }
            } 
            if(similarity==2){ //Double Dipper
              if(consim==1){//normal normal
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
                lgcont = gsimconNN(m0, v, s20, sumx, sumx2, xbar(p), njtmp, 1, 0, 1);
                lgconY = lgconY + lgcont;
              } 
              if(consim==2){
                //normal normal inverse gamma
                //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
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
            
            if(similarity == 1){ // Auxiliary
              //Dirichlet-Multinomial
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 0, 1);
              lgcatN = lgcatN + lgcatt;
            } 
            if(similarity==2){// Double dipper
              //Dirichlet-Multinomial
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 1, 1);
              lgcatN = lgcatN + lgcatt;
            } 
            
            njc(xcat(i*(ncat)+p)) += 1;
            njtmp += 1;
            
            if(similarity==1){
              //Dirichlet-Multinomial
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 0, 1);
              lgcatY = lgcatY + lgcatt;
            } 
            if(similarity==2){
              //Dirichlet-Multinomial
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcatt = gsimcatDM(njc, dirweights, catvec(p), 1, 1);
              lgcatY = lgcatY + lgcatt;
            } 
          }//chiude ciclo su p covariate discrete
          
          gtilY(j) = lgconY + lgcatY;
          gtilN(j) = lgconN + lgcatN;
          
          if(calibration == 1){
            
            ldSig = logdet(sigma_star_curr.row(j).t(), dim);
            weight(j) = dmvnorm(y.row(i).t(), mu_star_curr.row(j).t(),
                     sigma_star_curr.row(j).t(), dim, ldSig, 1) +
                       log((double) nj_curr(j)) + // cohesion part
                       lgcatY - lgcatN + // Categorical part
                       lgconY - lgconN;  // Continuous part
          } 
          if(calibration == 2){
            ldSig = logdet(sigma_star_curr.row(j).t(), dim);
            weight(j) = dmvnorm(y.row(i).t(), mu_star_curr.row(j).t(),
                     sigma_star_curr.row(j).t(), dim, ldSig, 1)+
                       log((double) nj_curr(j)) + // cohesion part
                       (1/((double)ncon + (double)ncat))*(lgcatY + lgconY - lgcatN - lgconN);
            }
          
      }

        //ASTERISCO
        // Uniform cohesion
        /*if(cohesion==2){
          ph(j) = ph(j) - log((double) nj_curr(j));
        }*/
      //}//chiude il ciclo sui j cluster
      for(j = nclu_curr; j < CC; j++){
        // Continuous Covariates
        lgcondraw = 0.0;
        //lgconN = 0.0;
        for(p = 0; p < (ncon); p++){
          tmp = xcon(i*(ncon) + p);
          
          if(similarity==1){ // Auxilliary
            if(consim==1){//normal normal
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcont = gsimconNN(m0, v, s20, tmp, tmp*tmp, xbar(p), 1, 0, 0, 1);
              lgcondraw += lgcont;
            } 
            if(consim==2){//normal normal inverse gamma
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcont = gsimconNNIG(m0, k0, nu0, s20, tmp, tmp*tmp, xbar(p), s2mle(p), 1, 0, 0, 1);
              lgcondraw += lgcont;
            } 
            
          } 
          if(similarity==2){ //Double Dipper
            if(consim==1){//normal normal
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
              lgcont = gsimconNN(m0, v, s20, tmp, tmp*tmp, xbar(p), 1, 1, 0, 1);
              lgcondraw += lgcont;
            } 
            if(consim==2){//normal normal inverse gamma
              //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
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
          
          if(similarity == 1){ // Auxiliary
            //Dirichlet-Multinomial
            //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
            lgcatt = gsimcatDM(njc, dirweights, catvec(p), 0, 1);
            lgcatdraw += lgcatt;
          } 
          if(similarity==2){// Double dipper
            //Dirichlet-Multinomial
            //vedi https://github.com/cran/ppmSuite/blob/master/src/Rutil.c
            lgcatt = gsimcatDM(njc, dirweights, catvec(p), 1, 1);
            lgcatdraw += lgcatt;
          }
        }//chiude ciclo su covariate discrete
        gtilY(j) = lgcondraw + lgcatdraw;
        gtilN(j) = lgcondraw + lgcatdraw;
        
        if(calibration == 1){
          ldSig = logdet(sigma_empty.row(j - nclu_curr).t(), dim);
          weight(j) = dmvnorm(y.row(i).t(), mu_empty.row(j - nclu_curr).t(),
                 sigma_empty.row(j - nclu_curr).t(), dim, ldSig, 1) +
                   log(alpha) - log(CC) +
                   lgcondraw + // Continuous covariate part
                   lgcatdraw; // categorical covariate part
          //lgcatY - lgcatN + // Categorical part
          //lgconY - lgconN;  // Continuous part;
          //}
        } 
        if(calibration == 2){
          ldSig = logdet(sigma_empty.row(j - nclu_curr).t(), dim);
          weight(j) = dmvnorm(y.row(i).t(), mu_empty.row(j - nclu_curr).t(),
                 sigma_empty.row(j - nclu_curr).t(), dim, ldSig, 1) +
                   log(alpha) - log(CC) +
                   (1/((double)ncon + (double)ncat))*(lgcondraw + lgcatdraw);
          //} 
        }
      }
      
      if(calibration == 1){
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
          
          ldSig = logdet(sigma_star_curr.row(j).t(), dim);
          weight(j) = dmvnorm(y.row(i).t(), mu_star_curr.row(j).t(),
                 sigma_star_curr.row(j).t(), dim, ldSig, 1) + 
                   log((double) nj_curr(j)) +  // Cohesion part
                   lgtilYk - lgtilNk; //This takes into account both cont and cat vars
          
          /*if(*cohesion==2){
           ph[k] = ph[k] - log((double) nh[k]);
          }*/
        }
        
        // calibration for a singleton
        for(j = nclu_curr; j < nclu_curr + CC; j++){
          ldSig = logdet(sigma_empty.row(j - nclu_curr).t(), dim);
          weight(j) = dmvnorm(y.row(i).t(), mu_empty.row(j - nclu_curr).t(),
                 sigma_empty.row(j - nclu_curr).t(), dim, ldSig, 1) +
                   log(alpha) - log(CC) +
                   lgtilN[j] - // Continuous covariate part
                   log(sgN);
        }
        
        /*if(*cohesion==2){// Note with a uniform cohesion, for a new cluster
         // the value of log(c({nclus_iter}}) = log(1) = 0;
         ph[nclus_iter] = ph[nclus_iter] - log(Mdp);
         }*/
        
      }
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
        weight(j) = weight(j)/denwei;
      }
      
      uu = R::runif(0.0,1.0);
      cweight = 0.0;
      
      for(j = 0; j < nclu_curr + CC; j++){
        cweight += weight(j);
        if (uu < cweight){
          newci = j;
          Rcpp::Rcout << "newci (986): " << newci << std::endl;
          break;
        }
      }
      if(newci<nclu_curr){//<=
        //Point 2a page 345 Favaro and Teh, second column
        curr_clu(i) = newci + 1;
        /*if(nj_curr(newci) == 0){
          nclu_curr += 1;
        }*/
        nj_curr(newci) += 1;
        Rcpp::Rcout << "RIGA 998ss " << std::endl;
        Rcpp::Rcout << "curr_clu: " << curr_clu << "nj_curr: " << nj_curr << std::endl;
        Rcpp::Rcout << "numero cluster: " << nclu_curr << std::endl;
      }else{
        id_empty = newci - nclu_curr;
        //### Point 2b  page 345 Favaro and Teh, second column
        curr_clu(i) = nclu_curr + 1;
        nclu_curr += 1;
        nj_curr(nclu_curr-1) = 1;///?
        //Rcpp::Rcout << "RIGA 986ss " << std::endl;
        //Rcpp::Rcout << "curr_clu: " << curr_clu << "nj_curr: " << nj_curr << std::endl;
        Rcpp::Rcout << "1008numerosità cluster: " << nj_curr << std::endl;
        /*lascio così perché C conta da 0
         * io voglio mettere 1 nella posizione nclu_curr+1
         * c(nj_current,1);*/
        //## the index of the empty component we have selected is
        //### update the phi parameters
        sigma_star_curr.row(nclu_curr-1) = sigma_empty.row(id_empty);
        mu_star_curr.row(nclu_curr-1) = mu_empty.row(id_empty);
        //## update empty parameter
        mu_empty.row(id_empty) = ran_mvnorm(hP0_m0, hP0_L0, dim).t();
        sigma_empty.row(id_empty) = ran_iwish(hP0_nu0, hP0_V0, dim).t();
      }
      }
      // Compute the CPO and lpml using the mixture
      //Rcpp::Rcout << "punto 000 - 999 " << std::endl;
      ldSig = logdet(sigma_star_curr.row(curr_clu(i)-1).t(), dim);
      like_iter(i) = dmvnorm(y.row(i).t(), mu_star_curr.row(curr_clu(i)-1).t(),
              sigma_star_curr.row(curr_clu(i)-1).t(), dim, ldSig, 1);
      if((l > (burn-1)) & (l % (thin) == 0)){
        
        // These are needed for WAIC
        mnlike(i) = mnlike(i) + (like_iter(i))/(double) nout;
        mnllike(i) = mnllike(i) + log(like_iter(i))/(double) nout;
        
        CPOinv(i) = CPOinv(i) + (1/(double) nout)*(1/like_iter(i));
      }
    }
    Rcpp::Rcout << "punto 00: " << l+1 << std::endl;
    //////////////////////////////////////////////////
    // update the cluster value with NEAL 8 (II step)
    //////////////////////////////////////////////////
    arma::vec Lnv(dim*dim);
    Rcpp::Rcout << "punto 0 " << std::endl;
    arma::mat Lnm(dim, dim);
    Rcpp::Rcout << "punto 1 " << std::endl;
    arma::mat L0mInv(dim, dim);
    Rcpp::Rcout << "punto 2 " << std::endl;
    arma::mat S0m(dim, dim);
    Rcpp::Rcout << "punto 3 " << std::endl;
    arma::mat Sigma(dim, dim);
    Rcpp::Rcout << "punto 4 " << std::endl;
    arma::mat nSigmaInv(dim, dim);
    Rcpp::Rcout << "punto 5 " << std::endl;
    arma::vec mun(dim);
    Rcpp::Rcout << "punto 6 " << std::endl;
    arma::vec theta(dim);
    arma::vec Snv(dim*dim);
    arma::mat Snm(dim, dim);
    arma::vec ybar(dim);
    ybar.fill(0.0);
    for(int row = 0; row < dim; row++){
      for(int col = 0; col < dim; col++){
        L0mInv(row, col) = hP0_L0(row * dim + col);
        S0m(row, col) = hP0_V0(row * dim + col);
      }
    }
    
    L0mInv = arma::inv(L0mInv);
    Rcpp::Rcout << "punto 1 " << std::endl;
    
    for(j = 0; j < nclu_curr; j++){
      arma::mat ytemp(nj_curr(j), dim);
      //Rcpp::Rcout << "RIGA 1041ss " << std::endl;
      //Rcpp::Rcout << "curr_clu: " << curr_clu << "nj_curr: " << nj_curr << std::endl;
      //Rcpp::Rcout << "numero cluster: " << nclu_curr << std::endl;
      int idxy = 0;
      for(ii = 0; ii < nobs; ii++){
        if(curr_clu(ii) == (j+1)){
          //Rcpp::Rcout << "curr_clu(ii)" << curr_clu(ii) << std::endl;
          ytemp.row(idxy) = y.row(ii);
          ybar += y.row(ii).t();
          idxy += 1;
        }
      }
      
      ybar = ybar/nj_curr(j);
      //Rcpp::Rcout << "ybar" << std::endl;
      Rcpp::Rcout << "punto 2 " << std::endl;
      for(int row = 0; row < dim; row++){
        for(int col = 0; col < dim; col++){
          
          Sigma(row, col) = sigma_star_curr(j, row * dim + col);
        }
      }
      Rcpp::Rcout << "some tedious stuff pt2" << std::endl;
      nSigmaInv = arma::inv(Sigma);
      nSigmaInv = nj_curr(j) * nSigmaInv;
      Lnm = arma::inv(L0mInv + nSigmaInv);//hP0_L0 + (nj_curr(j) * sigma_star_curr.row(j));
      mun = Lnm * (L0mInv * hP0_m0 + nSigmaInv * ybar);
      Rcpp::Rcout << "mun done" << std::endl;
      for(int row = 0; row < dim; row++){
        for(int col = 0; col < dim; col++){
          Lnv(row * dim + col) = Lnm(row, col);
          
        }
      }
    
     theta = ran_mvnorm(mun, Lnv, dim); 
      
      // Rcpp::Rcout << "theta done" << std::endl;
     Snm = arma::inv(S0m + ((ytemp.each_row()-theta.t()).t() * (ytemp.each_row()-theta.t())));
     
     for(int row = 0; row < dim; row++){
       for(int col = 0; col < dim; col++){
         Snv(row * dim + col) = Snm(row, col);
       }
     }
     // Rcpp::Rcout << "snm is hard" << std::endl;
     sigma_star_curr.row(j) = ran_iwish(hP0_nu0+nj_curr(j), Snv, dim).t();
     mu_star_curr.row(j) = ran_mvnorm(theta, sigma_star_curr.row(j).t(), dim).t();
     }

    //Rcpp::Rcout << "end of update" << std::endl;
    // in sample prediction to assess model fit
    //Rcpp::Rcout << "before in sample prediction" << std::endl;
    if((l > (burn-1)) & (l % (thin) == 0)){
      for(i = 0; i < nobs; i++){
        ispred_iter.row(i) = ran_mvnorm(mu_star_curr.row(curr_clu(i)-1).t(), 
                    sigma_star_curr.row(curr_clu(i)-1).t(), dim).t();
      }
    }
    //Rcpp::Rcout << "before saving" << std::endl;
    //////////////////////
    // Save MCMC iterates
    //////////////////////
    
    if((l > (burn-1)) & ((l + 1) % thin == 0)){
      
      nclus(ll) = nclu_curr;
      for(i = 0; i < nobs; i++){
        mu_out.row(ll*(nobs) + i) = mu_star_curr.row(curr_clu(i)-1);
        Rcpp::Rcout << "1" << std::endl;
        sigma_out.row(ll*(nobs) + i) = sigma_star_curr.row(curr_clu(i)-1);
        Rcpp::Rcout << "2" << std::endl;
        Clui(ll*(nobs) + i) = curr_clu(i);
        Rcpp::Rcout << "3" << std::endl;
        like(ll*(nobs) + i) = like_iter(i);
        Rcpp::Rcout << "4" << std::endl;
        //ispred.row(ll*(nobs) + i) = ispred_iter.row(i);
      }
      
      ll += 1;
    }
    Rcpp::Rcout << "saving" << std::endl;
    //Rcpp::Rcout << "curr_clu: " << curr_clu << "nj_curr: " << nj_curr << std::endl;
    //Rcpp::Rcout << "numero cluster: " << nclu_curr << std::endl;
  }//CLOSES MCMC iterations
  
  // calculate LPML
  /*lpml_iter=0.0;
  
  for(i = 0; i < nobs; i++){
    lpml_iter += log(1/CPOinv(i));
  }
  lpml = lpml_iter;
  */
  
  // Computing WAIC  (see Gelman article)
  
  /*elppdWAIC = 0.0;
  for(i = 0; i < nobs; i++){
    elppdWAIC += (2*mnllike(i) - log(mnlike(i)));
  }
  
  WAIC = -2*elppdWAIC;*/
  
  //RETURN
  return Rcpp::List::create(Rcpp::Named("mu") = mu_out,
                            Rcpp::Named("sigma") = sigma_out,
                            Rcpp::Named("nclu") = nclus,
                            Rcpp::Named("cl_lab") = Clui,
                            Rcpp::Named("like") = like,
                            Rcpp::Named("ispred") = ispred,
                            Rcpp::Named("nclus") = nclus);
                            //Rcpp::Named("WAIC") = WAIC,
                            //Rcpp::Named("lpml") = lpml
}
