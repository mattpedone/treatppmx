#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>

//#include <R_ext/Lapack.h>

namespace dens{
// Inverse Gamma density
// the parameterization here provides a mean of beta/(alpha - 1)
double dinvgamma(double y, double alpha, double beta, int logout){

  //	Rprintf("alpha = %f\n", alpha);
  //	Rprintf("beta = %f\n", beta);

  double ldens;

  ldens = alpha*log(beta) - lgamma(alpha) - (alpha + 1)*log(y) - (beta/y);

  if(logout) return ldens;
  return exp(ldens);

}

// Normal-Inverse Gamma density
// the parameterization here is such that mu0 is prior mean value with k0 prior "observations"
// and s20 is the prior guess for sig2 with nu0 prior "observations"
double dN_IG(double mu, double sig2, double mu0, double k0, double a0, double b0, int logout){

  //	Rprintf("alpha = %f\n", alpha);
  //	Rprintf("beta = %f\n", beta);

  double ldens;

  //	ldens =  0.5*(log(k0) - log(2*M_PI*sig2)) + a0*log(b0) -
  //	         lgammafn(a0) - (a0 + 1)*log(sig2) -
  //	         0.5*(1/sig2)*(k0*(mu-mu0)*(mu-mu0) + 2*b0);

  ldens = R::dnorm(mu, mu0, sqrt(sig2/k0),logout) +
    dens::dinvgamma(sig2, a0, b0, logout);
  //	Rprintf("ldens = %f\n", ldens);
  if(logout){ return ldens;
  }else{return exp(ldens);}

}
}//closes density namespace

namespace similarityf{
// normal-normal Similarity function with x following normal and m a normal as well (v is fixed).
// This is the result after integrating over mj.  One could also simply find the
// ratio between likelihood multiplied by prior divided by posterior

// The double dipper is included as an argument.

double gsimconNN(double m0, double v2, double s20, double sumx, double sumx2, double mle,
                 int n,  int DD, int cal, int logout){

  //	Rprintf("m0 = %f\n", m0);
  //	Rprintf("v2 = %f\n", v2);
  //	Rprintf("s20 = %f\n", s20);
  //	Rprintf("sumx = %f\n", sumx);
  //	Rprintf("sumx2 = %f\n", sumx2);
  //	Rprintf("n = %d\n", n);

  double mus, muss, s2s, s2ss;
  double ld1, ld2, ld3, ld4, ld5, ld6;
  double out;
  //	double out1, out2;

  s2s = 1/((n/v2) + (1/s20));
  mus = s2s*((1/v2)*sumx + (1/s20)*m0);

  //	Rprintf("mus = %f\n", mus);
  //	Rprintf("s2s = %f\n", s2s);

  s2ss = 1/((n/v2) + (1/s2s));
  muss = s2ss*((1/v2)*sumx + (1/s2s)*mus);

  //	Rprintf("muss = %f\n", muss);
  //	Rprintf("s2ss = %f\n", s2ss);

  ld1 = -0.5*n*log(2*M_PI*v2) - 0.5*(1/v2)*sumx2;
  ld2 = R::dnorm(m0, 0, sqrt(s20),1);
  ld3 = R::dnorm(mus, 0, sqrt(s2s),1);
  ld4 = R::dnorm(muss, 0, sqrt(s2ss),1);

  ld5 = R::dnorm(mle, m0, sqrt(s20),1);
  ld6 = R::dnorm(mle, mus, sqrt(s2s),1);
  //	Rprintf("ld1 = %f\n", ld1);
  //	Rprintf("ld2 = %f\n", ld2);
  //	Rprintf("ld3 = %f\n", ld3);
  //	Rprintf("ld4 = %f\n", ld4);


  out = ld1 + ld2 - ld3;
  //	out1 = ld1 + ld3 - ld4;
  //	out2 = ld5 - ld6;
  //	Rprintf("out = %f\n", out);
  //	Rprintf("out1 = %f\n", out1);
  if(DD==1) out = ld1 + ld3 - ld4;
  if(cal==1) out = ld5 - ld6;
  if(!logout) out = exp(out);
  return(out);

}

// normal-normal-IG Similarity function with x following normal and m,v a normal-IG.
// I didn't carry out integration explicitly over m and v. I simply used the fact that
// marginal likelihood (similarity) is equal to likelihood x prior / posterior.  This
// requires inserting a value for m and v which I use 0 and 1.

// The double dipper is included as an argument.


double gsimconNNIG(double m0, double k0, double nu0, double s20, double sumx, double sumx2,
                   double mnmle, double s2mle, int n, int DD, int cal, int logout){

  //	Rprintf("sumx = %f\n", sumx);
  //	Rprintf("sumx2 = %f\n", sumx2);
  //	Rprintf("n = %d\n", n);

  double a0, b0, m0s, m0ss, k0s, k0ss, a0s, a0ss, b0s, b0ss;
  double ld1, ld2, ld3, ld4, ld5, ld6, out;
  double mu=10, v2=0.1;
  double xbar = sumx*(1/ (double) n);

  a0 = 0.5*nu0;
  b0 = 0.5*nu0*s20;


  m0s = (k0*m0 + n*xbar)/(k0 + n);
  k0s = k0 + (double) n;
  a0s = a0 + 0.5*n;
  b0s = b0 + 0.5*(sumx2 - n*xbar*xbar) + 0.5*n*k0*(xbar - m0)*(xbar - m0)/(k0+n);

  m0ss = (k0s*m0s + n*xbar)/(k0s + n);
  k0ss = k0s + (double) n;
  a0ss = a0s + 0.5*n;
  b0ss = b0s + 0.5*(sumx2 - n*xbar*xbar) + 0.5*n*k0s*(xbar - m0s)*(xbar - m0s)/(k0s+n);

  ld1 = -0.5*n*log(2*M_PI*v2) - 0.5*(1/v2)*(sumx2 - 2*sumx*mu + n*mu*mu);
  ld2 = dens::dN_IG(mu, v2, m0, k0, a0, b0, 1);
  ld3 = dens::dN_IG(mu, v2, m0s, k0s, a0s, b0s, 1);
  ld4 = dens::dN_IG(mu, v2, m0ss, k0ss, a0ss, b0ss, 1);

  ld5 = dens::dN_IG(mnmle, s2mle, m0, k0, a0, b0, 1);
  ld6 = dens::dN_IG(mnmle, s2mle, m0s, k0s, a0s, b0s, 1);

  out = ld1 + ld2 - ld3;



  if(DD==1) out = ld1 + ld3 - ld4;
  if(cal==1) out = ld5 - ld6;
  if(!logout) out = exp(out);

  //	Rprintf("out = %f\n", out);
  return(out);

}

// Similarity function with for a categorical x  dirichlet-multinomial with out
// where only on object is allocated (x_i is basically univariate that identifies which
// category ith individual has).  The integral in reality is a product of two ratios,
// but one of the ratios is constant in terms of $x$ and so disappears in the ratio
// of similarity functions when updating cluster labels and so is ignored in the
// function that follows.
double gsimcatDM(arma::vec nobsj, arma::vec dirweights, int C, int DD, int logout){

  int ii, sumc;
  double tmp1=0.0,tmp2=0.0,tmp3=0.0,tmp4=0.0,tmp5=0.0,tmp6=0.0;
  double out;
  //	double ld1, ld2, ld3, ld4, out2, out3;
  //	double *pi_vec = R_Vector(C);
  //	double *dirweightsDD = R_Vector(C);
  //	double *dirweightsDDs = R_Vector(C);

  //	RprintIVecAsMat("nobsj", nobsj, 1, C);

  sumc=0;
  for(ii = 0; ii < C; ii++){
    sumc = sumc+nobsj(ii);

    tmp1 = tmp1 + dirweights(ii);
    tmp2 = tmp2 + lgamma(dirweights(ii));

    tmp3 = tmp3 + (double) nobsj(ii) + dirweights(ii);
    tmp4 = tmp4 + lgamma( (double) nobsj(ii) + dirweights(ii));

    tmp5 = tmp5 + 2*((double) nobsj(ii)) + dirweights(ii);
    tmp6 = tmp6 + lgamma( 2*((double) nobsj(ii)) + dirweights(ii));
  }

  //	Rprintf("tmp3 = %f\n", tmp3);
  //	Rprintf("tmp5 = %f\n", tmp5);

  // The next piece of code is for the similarity likelihoodxprior/posterior

  /*	for(ii=0;ii<C;ii++){
   sumc = sumc+nobsj[ii];
   pi_vec[ii] = 1/(double) C;
   dirweightsDD[ii] = (double) nobsj[ii] + dirweights[ii];
   dirweightsDDs[ii] = 2*((double) nobsj[ii]) + dirweights[ii];
  }
   ld1 = 0.0;
   for(ii=0;ii<C;ii++){
   ld1 = ld1 + nobsj[ii]*log(pi_vec[ii]);
   }
   ld2 = ddirich(pi_vec, dirweights, C, 1);
   ld3 = ddirich(pi_vec, dirweightsDD, C, 1);
   ld4 = ddirich(pi_vec, dirweightsDDs, C, 1);
   out1 = ld1 + ld2 - ld3;
   out2 = ld1 + ld3 - ld4;
   */

  out = (R::lgammafn(tmp1) - tmp2) + (tmp4 - R::lgammafn(tmp3));


  //	Rprintf("out = %f\n", out);

  if(DD==1) out = (R::lgammafn(tmp3) - tmp4) + (tmp6 - R::lgammafn(tmp5));
  if(sumc==0) out = log(1);
  if(!logout) out = exp(out);
  return(out);

}

}//closes namespace

// [[Rcpp::export]]
Rcpp::List myppmx(int iter, int burn, int thin, int nobs, int ncon, int ncat,
                  arma::vec catvec, double alpha, int CC, int cohesion,
                  int similarity, int consim, arma::vec y, arma::vec xcon,
                  arma::vec xcat, //int npred, arma::mat xconp, arma::mat xcatp,
                  arma::vec similparam, arma::vec modelpriors, arma::vec mhtune,
                  int calibration){

  //serve flag Reuse e m number of auxiliary parameters
  //fai check su eventuali conflitti nella fnuzione R

  // l - MCMC index
  // ll - MCMC index for saving iterates
  // i - individual index
  // ii - second individual index (for double for loops)
  // c - categorical variable index
  // p - number of covariates index
  // j - cluster index
  // t - subset of covariates index
  // mm - for auxiliary parameters
  int l, ll, i, ii, c, p, j, mm;

  double max_C, nout, sumx, sumx2;
  //number of saved iterations
  nout = (iter - burn)/(thin);
  Rcpp::Rcout << "nout =  " << nout << std::endl;
  //numero totale di covariate
  //int ncov = ncon + ncat;

  //maximum number of categories for categorical covariates
  max_C = catvec(0);

  for(c = 0; c < (ncat); c++){
    if(max_C < catvec[c]) max_C = catvec[c];
  }

  //sum of response vector
  double sumy = 0.0;
  for(i = 0; i < nobs; i++){
    sumy += y(i);
  }

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

  /////////////////////////////////
  // storing single MCMC iterate
  ////////////////////////////////

  double mu0_iter = 0.0;

  //upper bound on sig20
  //that is upper bound of uniform prior on cluster mean variance
  double sig20_iter = 0.5 * modelpriors(3);

  int nclus_iter = 0;
  int iaux;
  arma::vec curr_clu(nobs);
  curr_clu.fill(1);
  arma::vec nj_curr(nobs);
  nj_curr.fill(0);
  arma::vec njc((nobs)*(ncat));

  arma::vec xcontmp(nobs);

  //cfr adr1
  for(i = 0; i < nobs; i++){
    for(ii = 0; ii < nobs; ii++){
      if(curr_clu(i) == ii+1) nj_curr(ii) = nj_curr(ii) + 1;
    }
  }

  for(i = 0; i < nobs; i++){
    if(nj_curr(i) > 0) nclus_iter += 1;
  }


  arma::vec sig2h(nobs);
  sig2h.fill(0.5*modelpriors(3));
  arma::vec muh(nobs);
  muh.fill(0.0);

  /////////////////////////////////////
  // Storage for posterior predictive
  /////////////////////////////////////

  arma::vec ispred_iter(nobs);
  ispred_iter.fill(0.0);

  ////////////////////////////////////////
  // Storage needed to update parameters
  ////////////////////////////////////////

  // update Si (cluster labels);
  int njtmp;
  double auxm, auxs2, tmp;//, npdN, npdY, npd;
  double mudraw, sdraw, maxph, denph, cprobh, uu;
  double lgconN, lgconY, lgcatN, lgcatY, lgcondraw, lgcatdraw;
  double lgcont, lgcatt;

  arma::vec muaug(CC);
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


  // update muh (mustar)
  double mstar, s2star;

  // update mu0 (mean mustar)
  double summu;

  // update sig20 (variance of mustar)
  double llo, lln, llr, os0, ns0;

  // update sig2 (sigmastar)
  double osig, nsig;

  // Stuff for out of sample predictions
  //double lgcon0, lgcat0, mupred, sig2pred;

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

  //	Hyper-prior parameters

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

  for(mm = 0; mm < CC; mm++){
    muaug(mm) = R::rnorm(mu0_iter, sqrt(sig20_iter));
    saug(mm) = R::runif(smin, smax);
    }

  //storage for return
  arma::vec mu(nout * nobs, arma::fill::ones);
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

  ////////////////////////////////////////
  //
  // HERE COMES THE MCMC
  //
  ////////////////////////////////////////

  for(l = 0; l < iter; l++){

    if((l+1) % 10000 == 0){
      Rcpp::Rcout << "mcmc iter =  " << l+1 << std::endl;
    }

    //qui inizializza parametri P0 x reuse (vettori di lunghezza CC)

    /////////////////////////////////////////
    // update the cluster labels with NEAL 8
    /////////////////////////////////////////

    for(i = 0; i < nobs; i++){

      int zi = curr_clu(i)-1; //sottraggo 1 perché C conta da 0
      if(nj_curr(zi) > 1){
        /* l'osservazione NON è un singoletto
         * Neal lo chiama n_{-1, c}
         * nel modello di Raffaele S_{\tilde{j}}^{a\star, -i}
         */
        nj_curr(zi) = nj_curr(zi) - 1;
      }else{
        // l'osservazione è un singoletto
        iaux = curr_clu(i);
        if(iaux < nclus_iter){
          //questa condizione è vera se il singoletto "sta in mezzo" a non-singleton clusters
          //o meglio il singoletto non sta nell'ultimo cluster
          //quindi devo fare il relabel dei clusters
          //
          // faccio lo swap cluster labels tra curr_clu(i) e nclus_iter
          // ANCHE cluster specific parameters;

          // Metto tutti i membri dell'ultimo cluster nel cluster del soggetto i (iaux)
          for(ii = 0; ii < nobs; ii++){
            if(curr_clu(ii) == nclus_iter){
              curr_clu(ii) = iaux;
            }
          }

          //metto il soggetto i (singoletto) come ultimo cluster
          curr_clu(i) = nclus_iter;

          // The following steps swaps order of cluster specific parameters
          // so that the newly labeled subjects from previous step retain
          // their correct cluster specific parameters

          auxm = muh(iaux-1);
          muh(iaux-1) = muh(nclus_iter-1);
          muh(nclus_iter-1) = auxm;

          auxs2 = sig2h(iaux-1);
          sig2h(iaux-1) = sig2h(nclus_iter-1);
          sig2h(nclus_iter-1) = auxs2;


          // the number of members in cluster is also swapped with the last
          nj_curr(iaux-1) = nj_curr(nclus_iter-1);
          nj_curr(nclus_iter-1) = 1;
        }


        // Now remove the ith obs and last cluster;
        nj_curr(nclus_iter-1) = nj_curr(nclus_iter-1) - 1;
        nclus_iter = nclus_iter - 1;
      }

      // The atoms have been relabeled
      // Begin the cluster probabilities

      for(j = 0; j < nclus_iter; j++){

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
      for(mm = nclus_iter; mm < (nclus_iter+CC); mm++){
        gtilY(mm) = lgcondraw + lgcatdraw;
        gtilN(mm) = lgcondraw + lgcatdraw;
        }
      //gtilY(nclus_iter) = lgcondraw + lgcatdraw;
      //gtilN(nclus_iter) = lgcondraw + lgcatdraw;

      //EXT
      for(mm = nclus_iter; mm < (nclus_iter+CC); mm++){
        ph(mm) = R::dnorm(y(i), muaug(mm - nclus_iter), saug(mm - nclus_iter), 1) +
          log(alphadp) +  // DP part
          lgcondraw + // Continuous covariate part
          lgcatdraw; // categorical covariate part
        if(calibration == 2){
          ph(mm) = R::dnorm(y(i), muaug(mm - nclus_iter), saug(mm - nclus_iter), 1) +
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
          for(j = 1; j < nclus_iter + CC; j++){//1 - CC

          if(maxgtilN < gtilN(j)) maxgtilN = gtilN(j);

          if(j < nclus_iter){
            if(maxgtilY < gtilY(j)) maxgtilY = gtilY(j);
          }
        }

        sgY=0.0;
        sgN=0.0;

        for(j = 0; j < nclus_iter + CC; j++){//1 - CC

          lgtilN(j) = gtilN(j) - maxgtilN;
          sgN += exp(lgtilN(j));

          if(j < nclus_iter){// If x is included in an existing cluster in cannot be a singleton
            lgtilY(j) = gtilY(j) - maxgtilY;
            sgY = sgY + exp(lgtilY(j));
          }
        }

        //qui non modificare MA CONTROLLA!!
        // Calibrate the unnormalized cluster probabilities

        for(j = 0; j < nclus_iter; j++){
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
        for(mm = nclus_iter; mm < (nclus_iter+CC); mm++){
          ph(mm) = R::dnorm(y(i), muaug(mm - nclus_iter), saug(mm - nclus_iter), 1) +
            log(alphadp) +
            lgtilN(mm) - log(sgN) - log(CC);//in caso togli -log(CC)

          if(cohesion==2){// Note with a uniform cohesion, for a new cluster
            // the value of log(c({nclus_iter}}) = log(1) = 0;
            ph(mm) -= log(alphadp);
          }
        }

      }

      //NORMALIZZAZIONE PROBABILITà
      maxph = ph(0);
      for(j = 1; j < nclus_iter + CC; j++){//1 - CC
        if(maxph < ph(j)) maxph = ph(j);
      }

      denph = 0.0;
      for(j = 0; j < nclus_iter + CC; j++){//1 - CC
        ph(j) = exp(ph(j) - maxph);
        denph += ph(j);
      }

      for(j = 0; j < nclus_iter + CC; j++){//1 - CC
        probh(j) = ph(j)/denph;
      }

      uu = R::runif(0.0,1.0);

      //visto che al massimo posso aggiungere un solo cluster alla volta
      //lascerei +1 invece che + CC
      cprobh= 0.0;
      iaux = nclus_iter + 1;
      for(j = 0; j < nclus_iter + 1; j++){
        cprobh = cprobh + probh(j);
        if (uu < cprobh){
          iaux = j + 1;//arbitrariamente prendo il primo cfr Neal(2000)
          break;
        }
      }

      if(iaux <= nclus_iter){

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
        nclus_iter += 1;
        curr_clu(i) = nclus_iter;
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

    for(j = 0; j < nclus_iter; j++){
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
    for(j = 0; j < nclus_iter; j++){
      summu = summu + muh(j);
    }

    s2star = 1/(((double) nclus_iter/sig20_iter) + (1/s2));
    mstar = s2star*((1/sig20_iter)*summu + (1/s2)*m);

    mu0_iter = R::rnorm(mstar, sqrt(s2star));

    //sigma0 prior variance for muh
    os0 = sqrt(sig20_iter);
    ns0 = R::rnorm(os0,csigSIG0);

    if(ns0 > 0){

      lln = 0.0;
      llo = 0.0;
      for(j = 0; j < nclus_iter; j++){

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


    for(j = 0; j < nclus_iter; j++){
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

      nclus(ll) = nclus_iter;

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
