#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <stdlib.h>     /* srand, rand */
#include "utils_ct.h"

// [[Rcpp::export]]
Rcpp::List prior_ppmx_core(int iter, int burn, int thin, int nobs,
                          int PPMx, int ncon, int ncat, double alpha,
                          double sigma, arma::mat Vwm, int cohesion, int CC,
                          int consim, int similarity,
                          int calibration, int coardegree, arma::vec xcon,
                          arma::vec similparam,
                          arma::vec curr_cluster,
                          arma::vec card_cluster, int ncluster_curr){

  int l, ll, i, ii, p, j, jj, zi;

  int nout = (iter - burn)/(thin); //number of saved iterations

  //////////////////////////
  ////Cluster-related stuff
  //////////////////////////

  arma::vec curr_clu = curr_cluster;
  arma::vec nj_curr(nobs);
  for(i = 0; i < nobs; i++){
    nj_curr(i) = card_cluster(i);
  }

  int nclu_curr = ncluster_curr;
  arma::vec nclu(nout, arma::fill::zeros);
  arma::mat njout(nobs, nout, arma::fill::zeros);

  arma::vec weight(nobs + CC);
  arma::vec pweight(nobs + CC);
  double vp;
  vp = pow(nobs + CC, -1.0);
  weight.fill(vp);

  ////////////////////////////////////////////////////
  //// cluster specific suffiient statistics
  ////////////////////////////////////////////////////

  arma::vec sumx(nobs * ncon, arma::fill::zeros);
  arma::vec sumx2(nobs * ncon, arma::fill::zeros);
  double auxreal, sumxtmp, sumx2tmp;

  //arma::vec auxv(dim);
  int iaux;
  if(PPMx == 1){
    // Fill in cluster-specific sufficient statistics based on first partition
    for(i = 0; i < nobs; i++){
      for(p = 0; p < ncon; p++){
        sumx((curr_clu(i)-1)*ncon + p) += xcon(i * ncon + p);
        sumx2((curr_clu(i)-1)*ncon + p) += xcon(i * ncon + p) * xcon(i * ncon + p);
      }
    }
  }

  ////////////////////////////////////////
  // Stuff needed for similarities
  ////////////////////////////////////////

  double lgconN, lgconY, lgcatN, lgcatY, tmp;
  double lgcont;

  double m0 = similparam(0);
  double s20 = similparam(1);
  double v = similparam(2);
  double k0 = similparam(3);
  double nu0 = similparam(4);

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
  // Stuff needed for probabilities
  ////////////////////////////////////////
  double maxwei, denwei, uu, cweight;
  int newci, id_empty;

  ll = 0;

  ////////////////////////////////////////
  //
  // HERE COMES THE MCMC
  //
  ////////////////////////////////////////

  for(l = 0; l < iter; l++){


    //it = 0;
    for(i = 0; i < nobs; i++){


      /////////////////////////////////////////
      // update the cluster labels with NEAL 8
      /////////////////////////////////////////
      zi = curr_clu(i)-1; //sottraggo 1 perché C conta da 0

      if(nj_curr(zi) > 1){// Case where the nj corresponding to zi is such that nj>1
        nj_curr(zi) -= 1;
        if(PPMx == 1){
          // need to reduce sumx sumx2
          for(p = 0; p < ncon; p++){
            sumx(zi * ncon + p) -= xcon(i * ncon + p);
            sumx2(zi * ncon + p) -= xcon(i * ncon + p) * xcon(i * ncon + p);
          }
        }
      } else {// Case where the nj corresponding to zi is such that nj=1

        iaux = curr_clu(i);
        if(iaux < nclu_curr){
          for(ii = 0; ii < nobs; ii++){
            if(curr_clu(ii) == nclu_curr){
              curr_clu(ii) = iaux;
            }
          }
          curr_clu(i) = nclu_curr;

          nj_curr(iaux - 1) = nj_curr(nclu_curr - 1);
          nj_curr(nclu_curr - 1) = 1;

          if(PPMx == 1){
            // need to swap sumx and sumx2
            for(p = 0; p < ncon; p++){
              auxreal = sumx((iaux - 1)*(ncon) + p);
              sumx((iaux-1)*ncon + p) = sumx((nclu_curr-1)*(ncon) + p);
              sumx((nclu_curr-1)*(ncon) + p) = auxreal;

              auxreal = sumx2((iaux-1)*(ncon) + p);
              sumx2((iaux-1)*(ncon) + p) = sumx2((nclu_curr-1)*(ncon) + p);
              sumx2((nclu_curr-1)*(ncon) + p) = auxreal;
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
        }

        nclu_curr -= 1;


        //FINE PARTE DELICATA
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

            if(similarity==1){ // Auxilliary
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, nj_curr(j), 0, 1);
                lgconN += lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, nj_curr(j), 0, 1);
                lgconN += lgcont;
              }
            }
            if(similarity==2){ //Double Dipper
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, nj_curr(j), 1, 1);
                lgconN += lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, nj_curr(j), 1, 1);
                lgconN += lgcont;
              }
            }

            // now add ith individual back;
            sumxtmp += xcon(i*(ncon)+p);
            sumx2tmp += xcon(i*(ncon)+p)*xcon(i*(ncon)+p);

            if(similarity==1){ // Auxilliary
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, nj_curr(j) + 1, 0, 1);
                lgconY += lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, nj_curr(j) + 1, 0, 1);
                lgconY += lgcont;
              }
            }
            if(similarity==2){ //Double Dipper
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, nj_curr(j) + 1, 1, 1);
                lgconY += lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, nj_curr(j) + 1, 1, 1);
                lgconY += lgcont;
              }
            }
          }//chiude ciclo su p covariate continue

          gtilY(j) = lgconY + lgcatY;
          gtilN(j) = lgconN + lgcatN;


        }// this closes PPMx

        if((PPMx == 0) | ((calibration != 2) & (PPMx == 1))){
          // cohesion part
          if(cohesion == 1){
            weight(j) = log((double) nj_curr(j));
          }
          if(cohesion == 2){
            //weight(j) = log((double) (Vwm(nobs-1, nclu_curr-1)/Vwm(nobs-2, nclu_curr-1))+(nj_curr(j)-sigma));
            weight(j) = log((double) Vwm(nobs-1, nclu_curr-1)) - log((double) Vwm(nobs-2, nclu_curr-1));
            weight(j) += lgamma(nj_curr(j) + 1 - sigma) - lgamma(nj_curr(j) - sigma);

            //Rcpp::Rcout << "weight: " << weight(tt, j) << std::endl;
          }
          weight(j) += lgcatY - lgcatN + // Categorical part
            lgconY - lgconN;  // Continuous;
        }

        if((calibration == 2) & (PPMx == 1)){
          //ast: qui ci metto l'if per coehesion
          // cohesion part
          if(cohesion == 1){
            weight(j) = log((double) nj_curr(j));
          }
          if(cohesion == 2){
            //weight(j) = log((double) (Vwm(nobs-1, nclu_curr-1)/Vwm(nobs-2, nclu_curr-1))+(nj_curr(j)-sigma));
            weight(j) = log((double) Vwm(nobs-1, nclu_curr-1)) - log((double) Vwm(nobs-2, nclu_curr-1));
            weight(j) += lgamma(nj_curr(j) + 1 - sigma) - lgamma(nj_curr(j) - sigma);
          }
          if(coardegree == 1){
            weight(j) += (1/((double)ncon + (double)ncat))*(lgcatY + lgconY - lgcatN - lgconN);
          }
          if(coardegree == 2){
            weight(j) += (1/(pow(((double)ncon + (double)ncat), 1.0/2.0)))*(lgcatY + lgconY - lgcatN - lgconN);
          }
        }
      } //chiude la similarity & calibration for current clusters


      for(j = nclu_curr; j < (nclu_curr + CC); j++){
        jj = j - nclu_curr;
        lgcondraw = 0.0;
        lgcatdraw = 0.0;
        if(PPMx == 1){
          // Continuous Covariates
          for(p = 0; p < (ncon); p++){
            tmp = xcon(i*(ncon) + p);
            if(similarity==1){ // Auxilliary
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, tmp, tmp*tmp, 1, 0, 1);
                lgcondraw += lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, tmp, tmp*tmp, 1, 0, 1);
                lgcondraw += lgcont;
              }
            }
            if(similarity==2){ //Double Dipper
              if(consim==1){//normal normal
                lgcont = gsimconNN(m0, v, s20, tmp, tmp*tmp, 1, 1, 1);
                lgcondraw += lgcont;
              }
              if(consim==2){//normal normal inverse gamma
                lgcont = gsimconNNIG(m0, k0, nu0, s20, tmp, tmp*tmp, 1, 1, 1);
                lgcondraw += lgcont;
              }
            }
          }
          gtilY(j) = lgcondraw + lgcatdraw;
          gtilN(j) = lgcondraw + lgcatdraw;//ATTENZIONE
        }//closes PPMx

        if((PPMx == 0) | ((calibration != 2) & (PPMx == 1))){
          // cohesion part
          if(cohesion == 1){
            weight(j) = log(alpha) - log(CC);
          }
          if(cohesion == 2){
            //weight(j) = log((double) (Vwm(nobs-1, nclu_curr)/Vwm(nobs-2, nclu_curr-1)));//*(nj_curr(tt, j)-sigma);
            weight(j) = log((double) Vwm(nobs-1, nclu_curr-1)) - log((double) Vwm(nobs-2, nclu_curr-1));
            //weight(j) += lgamma(nj_curr(j) + 1 - sigma) - lgamma(nj_curr(j) - sigma);
          }
          weight(j) += lgcondraw + // Continuous covariate part
            lgcatdraw; // categorical covariate part
        }

        if((calibration == 2) & (PPMx == 1)){
          // cohesion part
          if(cohesion == 1){
            weight(j) = log(alpha) - log(CC);
          }
          if(cohesion == 2){
            //weight(j) = log((double) (Vwm(nobs-1, nclu_curr)/Vwm(nobs-2, nclu_curr-1)));//*(nj_curr(tt, j)-sigma);
            weight(j) = log((double) Vwm(nobs-1, nclu_curr-1)) - log((double) Vwm(nobs-2, nclu_curr-1));
            //weight(j) += lgamma(nj_curr(j) + 1 - sigma) - lgamma(nj_curr(j) - sigma);
          }
          if(coardegree == 1){
            weight(j) += (1/((double)ncon + (double)ncat))*(lgcondraw + lgcatdraw);
          }
          if(coardegree == 2){
            weight(j) += (1/(pow(((double)ncon + (double)ncat), 1.0/2.0)))*(lgcondraw + lgcatdraw);
          }
        }
      }//this ends loop on empty clusters

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
            sgY += exp(lgtilY(j));
          }
        }

        // Calibrazione prob di cluster esistenti
        for(j = 0; j < nclu_curr; j++){
          lgtilNk = lgtilN(j) - log(sgN);
          lgtilYk = lgtilY(j) - log(sgY);

          //ast: qui ci metto if per cohesion (dentro perché dipende da j)
          // cohesion part
          if(cohesion == 1){
            weight(j) = log((double) nj_curr(j));
          }
          if(cohesion == 2){
            //weight(j) = log((double) (Vwm(nobs-1, nclu_curr-1)/Vwm(nobs-2, nclu_curr-1))+(nj_curr(j)-sigma));
            weight(j) = log((double) Vwm(nobs-1, nclu_curr-1)) - log((double) Vwm(nobs-2, nclu_curr-1));
            weight(j) += lgamma(nj_curr(j) + 1 - sigma) - lgamma(nj_curr(j) - sigma);
          }
          weight(j) += lgtilYk - lgtilNk; //cov cont and cat
        }

        // calibration for empty clusters
        for(j = nclu_curr; j < nclu_curr + CC; j++){
          jj = j - nclu_curr;
          lgtilNk = lgtilN(j) - log(sgN);

          // cohesion part
          if(cohesion == 1){
            weight(j) = log(alpha) - log(CC);
          }
          if(cohesion == 2){
            //weight(j) = log((double) (Vwm(nobs-1, nclu_curr)/Vwm(nobs-2, nclu_curr-1)));//*(nj_curr(tt, j)-sigma);
            weight(j) = log((double) Vwm(nobs-1, nclu_curr-1)) - log((double) Vwm(nobs-2, nclu_curr-1));
            //weight(j) += lgamma(nj_curr(j) + 1 - sigma) - lgamma(nj_curr(j) - sigma);
          }
          weight(j) += lgtilNk;
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
        //Rcpp::Rcout << "peso cluster " << j << ": " << weight(j) << std::endl;
        denwei += weight(j);
      }
      for(j = 0; j < nclu_curr + CC; j++){
        pweight(j) = weight(j)/denwei;
        //Rcpp::Rcout << "peso cluster " << j << ": " << pweight(j) << std::endl;
      }

      //sample the new cluster for i-th observation
      uu = R::runif(0.0,1.0);
      cweight = 0.0;
      newci = nclu_curr + CC;
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
        nj_curr(curr_clu(i) - 1) += 1;

        /*for(k = 0; k < dim; k++){
         loggamma(i, k) = calculate_gamma(eta_star_curr.slice(tt),
         z, beta, curr_clu(tt, it)-1, k, i, 1);
        }*/
      }else{
        id_empty = newci - nclu_curr - 1;
        nclu_curr += 1;
        curr_clu(i) = nclu_curr;
        nj_curr(nclu_curr-1) = 1;
      }
      if(PPMx == 1){
        // need to now add the xcon to the cluster to which it was assigned;
        for(p = 0; p < ncon; p++){
          sumx((curr_clu(i)-1)*ncon + p) += xcon(i * ncon + p);
          sumx2((curr_clu(i)-1)*ncon + p) += xcon(i * ncon + p) * xcon(i * ncon + p);
        }
      }
    }//this closes the loop on observations i

    //this closes the loop on candidate therapies T
    if((l > (burn-1)) & ((l + 1) % thin == 0)){
      nclu(ll) = nclu_curr;
      njout.col(ll) = nj_curr;
      ll += 1;
    }
  }//CLOSES MCMC iterations

  //RETURN
  return Rcpp::List::create(//Rcpp::Named("mu") = mu_out,
    Rcpp::Named("nclu") = nclu,
    Rcpp::Named("njout") = njout);
}
