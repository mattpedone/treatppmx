#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>

/*
 * functions quform, inner_product and squared norm are adapted from
 * correspondent functions in MATRIX.C - Routines for doing matrix operations.
 *
 *  Copyright (c) 1995-2004 by Radford M. Neal
 *
 * Permission is granted for anyone to copy, use, modify, or distribute this
 * program and accompanying programs and documents for any purpose, provided
 * this copyright notice is retained and prominently displayed, along with
 * a note saying that the original programs are available from Radford Neal's
 * web page, and note is made of any changes made to the programs.  The
 * programs and documents are distributed without any warranty, express or
 * implied.  As the programs were written for research purposes only, they have
 * not been tested to the degree that would be advisable in any important
 * application.  All use of these programs is entirely at the user's own risk.
 */

/* most of remaining are adapted from correspondent function in ppmSuite package
 * Copyright (c) 2009 Steven McKay Curtis and Garritt L. Page
 *
 * I give permission for anyone to use and/or modify these
 * programs provided the following conditions are met:
 *
 * 1) This copyright message is retained.
 * 2) Any changes to the code are noted.
 *
 * These programs are provided WITHOUT ANY WARRNTY.
 *
 */

//quadratic form
double quform(arma::vec x, arma::vec A, int dim){

  int i, j;
  double sm = 0.0;
  for(i = 1; i < dim; i++){
    for(j = 0; j < i; j++){
      sm += x(i)*x(j)*A(i * dim + j);
    }
  }
  sm *= 2;
  for(i = 0; i < dim; i++){
    sm += x(i)*x(i)*A(i * dim + i);
  }
  return sm;
}

//inner product
/* Each vector is of length n, and is stored in memory in successive locations
 * that are at a given distance apart. For an ordinary vector, a distance of 1
 * is used, but other distances can come about from looking at columns of a
 * matrix.
 */

double inner_product(arma::vec v1, int d1, arma::vec v2, int d2, int n){
  double s = 0;
  int i;
  int idx1= 0;
  int idx2 = 0;
  for(i = n; i > 0; i--){
    s += v1(idx1) * v2(idx2);
    idx1 += d1;
    idx2 += d2;
  }
  return s;
}

//squared norm
/* The vector (of length n) is stored in memory in successive locations that are
 * at a given distance apart.  For an ordinary vector, a distance of 1 is
 * used, but other distances can come about from looking at columns of a matrix.
 */

double squared_norm(arma::vec v, int d, int n, int sq){
  double s = 0;
  int i;
  int idx = 0;

  for(i = n; i > 0; i--){
    s += v(idx)*v(idx);
    idx += d;
  }
  if(sq == 0){
    s = std::pow(s, .5);
    return s;
  }
  return s;
}

//Cholesky decomposition of posdef sym n x n matrix
arma::vec cholesky(arma::vec A, int n) {
  arma::vec L(n*n);
  L.fill(0);

  for (int i = 0; i < n; i++){
    for(int j = 0; j < (i+1); j++){
      double s = 0;
      for(int k = 0; k < j; k++){
        s += L[i * n + k] * L[j * n + k];
      }
      if(i == j){
        L[i * n + j] = sqrt(A[i * n + i] - s);
      } else {
        L[i * n + j] = (1.0 / L[j * n + j] * (A[i * n + j] - s));
      }
    }
  }
  return L;
}

// [[Rcpp::export]]
double myround( double x )
{
  const double sd = 1000; //for accuracy to 3 decimal places
  return int(x*sd + (x<0? -0.5 : 0.5))/sd;
}

/*
 * If the symmetric positive definite matrix A is represented by its Cholesky
 * decomposition A = LL' or A = U'U, then the determinant of this matrix can
 * be calculated as the product of squares of the diagonal elements of L or U.
 * vector A of dim (n x n) x 1 stores the cholesky decomposition of the
 * (n x n) matrix whose log determinant we are looking for
 * A n x n matrix whose determinant we want to compute
 */
double logdet(arma::vec A, int n) {
  arma::vec L(n*n);
  L.fill(0);
  double logdet = 0.0;
  double de = 1.0;

  for (int i = 0; i < n; i++){
    for(int j = 0; j < (i+1); j++){
      double s = 0;
      for(int k = 0; k < j; k++){
        s += L[i * n + k] * L[j * n + k];
      }
      if(i == j){
        L[i * n + j] = sqrt(A[i * n + i] - s);
      } else {
        L[i * n + j] = (1.0 / L[j * n + j] * (A[i * n + j] - s));
      }
      if(i == j){
        de *= (A[i * n + i] - s);
      }
      logdet =log(de);
    }
  }
  return logdet;
}

// Inverse Gamma density
// the parameterization here provides a mean of beta/(alpha - 1)
double dinvgamma(double y, double alpha, double beta, int logout){

  //	Rprintf("alpha = %f\n", alpha);
  //	Rprintf("beta = %f\n", beta);

  double ldens, out;

  ldens = alpha*log(beta) - lgamma(alpha) - (alpha + 1)*log(y) - (beta/y);

  if(logout){
    out = ldens;
  } else{
    out = exp(ldens);
  }
  return out;
}

// Normal-Inverse Gamma density
// the parameterization here is such that mu0 is prior mean value with k0 prior "observations"
// and s20 is the prior guess for sig2 with nu0 prior "observations"
double dN_IG(double mu, double sig2, double mu0, double k0, double a0, double b0, int logout){

  //	Rprintf("alpha = %f\n", alpha);
  //	Rprintf("beta = %f\n", beta);

  double ldens, out;

  //	ldens =  0.5*(log(k0) - log(2*M_PI*sig2)) + a0*log(b0) -
  //	         lgammafn(a0) - (a0 + 1)*log(sig2) -
  //	         0.5*(1/sig2)*(k0*(mu-mu0)*(mu-mu0) + 2*b0);

  ldens = R::dnorm(mu, mu0, sqrt(sig2/k0),logout) + dinvgamma(sig2, a0, b0, logout);
  //	Rprintf("ldens = %f\n", ldens);
  if(logout){
    out = ldens;
  }else{
    out = exp(ldens);
  }
  return out;
}

/* Multivariate normal density
 * y is the observation vector for which the density will be computed
 * mu is mean vector
 * Sig è la VARIANCE matrix
 * dim dimension of multivariate distribution
 * ld vuole log determinant variance matrix
 * logout logical if 1 returns log density
 */
double dmvnorm(arma::vec y, arma::vec mu, arma::vec Sig, int dim, double ld,
               int logout){

  int i, j, k;
  arma::mat work(dim, dim);

  for(j = 0; j < dim; j++){
    for(k = 0; k < dim; k++){
      work(j, k) = Sig(j * dim + k);
    }
  }
  work = arma::inv(work);

  for(j = 0; j < dim; j++){
    for(k = 0; k < dim; k++){
      Sig(j * dim + k) = work(j, k);
    }
  }

  arma::vec scr(dim);

  double qf, out;
  for(i = 0; i < dim; i++){
    scr(i) = y(i) - mu(i);
  }
  qf = quform(scr, Sig, dim);
  out = -(double) dim*M_LN_SQRT_2PI - 0.5*(ld + qf);
  if (logout)	return out;
  return exp(out);
}

/* Sampling Multivariate Normal
 * m vector of dimension dim storing the mean
 * Sig vector for VARIANCE MATRIX
 * dim dimension of multivarite normal distribution
 */
arma::vec ran_mvnorm(arma::vec m, arma::vec Sig, int dim){
  int i,j;
  arma::vec cholV(dim * dim);
  cholV = cholesky(Sig, dim);
  arma::vec z(dim);
  arma::vec out(dim);
  for(i = 0; i < dim; i++){
    z(i) = R::rnorm(0,1);
    out(i) = m(i);
    for(j = 0; j <= i; j++){
      out(i) += cholV(i * dim + j)*z(j);
    }
  }
  return out;
}

/* The following provides a density function for the Inverse Wishart
 * distribution.
 * Sig is argument and S is parameter.
 * SSiginv - is SSig^{-1} that is found in the trace of the density stored in a vec
 * In Hoff is S_0\Sigma^{-1}
 * dim id the dimension of the scale matrix distribution
 * detSig is the determinant of Sig (argument)
 * detS is the determinant of S (parameter)
 * nu0 - degrees of freedom of the inverse-wishart function
 * logout logical. if 1 returns logdensity
 */
double dinvwish(arma::vec SSiginv, int dim, double detSig, double detS, int nu0,
                int logout){

  int i;
  double out;
  double lgammasum = 0.0;
  double trace = 0.0;
  double p1, p2, p3, p4;

  for(i = 0; i < dim; i++){
    lgammasum += R::lgammafn(0.5*(nu0 + 1 - (i+1))); //I have to add 1 since C starts at zero
  }

  for(i = 0; i < dim*dim; i++){
    if(i % (dim+1) == 0){
      trace += SSiginv(i);
    }
  }
  p1 = 0.5 * nu0 * dim;
  p2 = 0.25 * dim * (dim - 1);//https://it.wikipedia.org/wiki/Coefficiente_binomiale#Caso_particolare
  p3 = 0.5 * nu0;
  p4 = 0.5 * (nu0 + dim + 1);

  out = -p1*log(2) - (p2*log(M_PI) + lgammasum) + p3*log(detS) - p4*log(detSig) - 0.5*trace;
  if(logout) return out;

  return exp(out);
}

/* Random Draw from Wishart distribution
 * nu degrees of freedom
 * Sig matrix parameter NOT THE CHOLESKY DECOMPOSITION
 * dim dimension of the Wishart distribution
 */
arma::vec ran_iwish(int nu, arma::vec Sig, int dim){

  int i, j, k;

  arma::vec x(dim);
  arma::vec zeros(dim);
  zeros.fill(0.0);
  //arma::vec outv(dim * dim);
  arma::vec out(dim * dim);
  out.fill(0.0);

  for(i = 0; i < nu; i++){
    x = ran_mvnorm(zeros, Sig, dim);
    for(j = 0; j < dim; j++){
      for(k = 0; k <= j; k++){
        out(j * dim + k) += x(j)*x(k);
      }
    }
  }

  /* fill the upper triangular part with lower triangular part */
  for(j=0; j<dim; j++){
    for(k=0; k<j; k++){
      out(k * dim + j) = out(j * dim + k);
    }
  }

  arma::mat outmat(dim, dim);
  for(j = 0; j < dim; j++){
    for(k = 0; k < dim; k++){
      outmat(j,  k) = out(j * dim + k);
    }
  }
  outmat = arma::inv(outmat);

  for(j = 0; j < dim; j++){
    for(k = 0; k < dim; k++){
      out(j * dim + k) = outmat(j,  k);
    }
  }
  return out;
}

/* normal-normal Similarity function with x following normal and m a normal as well (v is fixed).
 * This is the result after integrating over mj.  One could also simply find the
 * ratio between likelihood multiplied by prior divided by posterior
 * The double dipper is included as an argument.
 */
double gsimconNN(double m0, double v2, double s20, double sumx, double sumx2,
                 int n,  int DD, int logout){

  double mus, muss, s2s, s2ss;
  double ld1, ld2, ld3, ld4, ld5, ld6;
  double out;

  s2s = 1/((n/v2) + (1/s20));
  mus = s2s*((1/v2)*sumx + (1/s20)*m0);

  s2ss = 1/((n/v2) + (1/s2s));
  muss = s2ss*((1/v2)*sumx + (1/s2s)*mus);

  ld1 = -0.5*n*log(2*M_PI*v2) - 0.5*(1/v2)*sumx2;
  ld2 = R::dnorm(m0, 0, sqrt(s20),1);
  ld3 = R::dnorm(mus, 0, sqrt(s2s),1);
  ld4 = R::dnorm(muss, 0, sqrt(s2ss),1);

  out = ld1 + ld2 - ld3;
  if(DD==1) out = ld1 + ld3 - ld4;

  if(!logout) out = exp(out);
  return out;
}

/* normal-normal-IG Similarity function with x following normal and m,v a normal-IG.
 * I didn't carry out integration explicitly over m and v. I simply used the fact that
 * marginal likelihood (similarity) is equal to likelihood x prior / posterior.  This
 * requires inserting a value for m and v which I use 0 and 1.
 * The double dipper is included as an argument.
 */
double gsimconNNIG(double m0, double k0, double nu0, double s20, double sumx, double sumx2,
                   int n, int DD, int logout){

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
  ld2 = dN_IG(mu, v2, m0, k0, a0, b0, 1);
  ld3 = dN_IG(mu, v2, m0s, k0s, a0s, b0s, 1);
  ld4 = dN_IG(mu, v2, m0ss, k0ss, a0ss, b0ss, 1);

  out = ld1 + ld2 - ld3;

  if(DD==1) out = ld1 + ld3 - ld4;
  if(!logout) out = exp(out);

  return out;
}

/* Similarity function with for a categorical x  dirichlet-multinomial with out
 * where only on object is allocated (x_i is basically univariate that identifies which
 * category ith individual has).  The integral in reality is a product of two ratios,
 * but one of the ratios is constant in terms of $x$ and so disappears in the ratio
 * of similarity functions when updating cluster labels and so is ignored in the
 * function that follows.
 */
double gsimcatDM(arma::vec nobsj, arma::vec dirweights, int C, int DD, int logout){

  int ii, sumc;
  double tmp1=0.0,tmp2=0.0,tmp3=0.0,tmp4=0.0,tmp5=0.0,tmp6=0.0;
  double out;

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

  out = (R::lgammafn(tmp1) - tmp2) + (tmp4 - R::lgammafn(tmp3));

  if(DD==1) out = (R::lgammafn(tmp3) - tmp4) + (tmp6 - R::lgammafn(tmp5));
  if(sumc==0) out = log(1);
  if(!logout) out = exp(out);

  return out;
}

double calculate_gamma(arma::mat eta, arma::mat ZZ, arma::vec beta, int clu_lg,
                       int k, int i, int Log){
  /* function that compute log linear predictor.
   * inputs:
   * - eta: mvn intercept (BNP)
   * - ZZ: prognostic covariates
   * - beta: prognostic coefficients
   * - clu_lg: cluster of loggamma we are currently working on
   * - k: index for category
   * - i: index for subject
   * - Log if == 1 returns log gamma
   */

  int q, h; //indices for prognostic covariates (real &"instrumental")
  int Q = ZZ.n_cols;
  int dim = eta.n_rows;

  double lg = 0.0;
  double gamma_ik;

  if(eta.n_cols == 1){
    lg = eta(k);
  } else {
    lg = eta(k, clu_lg);
  }

  for(q = 0; q < Q; q++){
    h = k + q * dim;
    //h = q + k * Q;
    lg += (beta(h) * ZZ(i, q));
  }

  if(Log == 1){
    gamma_ik = lg;
  } else{
    gamma_ik = exp(lg);
  }
  return gamma_ik;
}

arma::vec ran_dir(arma::vec param, double delta, double rho){
  int J = param.size();
  arma::vec dir_sample(J);
  dir_sample.zeros();
  double sum = 0;
  arma::vec param_scale_shift(J);
  param_scale_shift = param*rho + delta;

  for( int j = 0; j < J; j++){
    dir_sample[ j ] = R::rgamma(param_scale_shift(j), 1);//[0];
    sum += dir_sample[j];
  }

  // Protect against dividing by zero
  if(sum == 0){
    sum = 1;
  }

  dir_sample = dir_sample/sum;

  return dir_sample;
}

double log_mult(arma::mat y, arma::mat JJ){
  //multinomial probability mass function expressed using gamma function
  int nobs, i;
  double log_mult = 0.0;
  nobs = JJ.n_rows;

  for(i = 0; i < nobs; i++){
    log_mult += lgamma(arma::sum(y.row(i)) + 1) +
      arma::sum(-lgamma(y.row(i) + 1) + y.row(i)%log(JJ.row(i)));
  }

  return log_mult;
}

/*
 * the following function updates the mvn intercept clusterwise
 * $\boldsymbol{\eta}_{j}^{\star}$
 */

Rcpp::List eta_update(arma::mat JJ, arma::mat loggamma,
                      arma::vec curr_clu, //arma::vec nj_curr,
                      arma::vec treatments, int t,
                      arma::vec eta, arma::vec eta_flag,
                      arma::vec mu_star, arma::vec sigma_star, int jj, double mhtune){
  /*
   * JJ: matrix of independent gamma variables
   * loggamma: matrix of log-linear predictor
   * curr_clu: vectors of labels of current cluster assignment
   * nj_curr: vectors of number of individuals in each cluster
   * eta: vector of mvn intercept of the jj-th cluster
   * eta_flag: flag for eta MH acceptance
   * mu_star: vector of jj-th cluster prior mean
   * sigma_star: vector of jj-th cluster prior covariance (matrix stored in array)
   * jj: cluster we are working on
   */

  int dim = JJ.n_cols;
  int nobs = JJ.n_rows;

  int i, k, it;//, ii;
  /*
   * indices for:
   * i: individuals
   * k: categories
   * ii: second index for individuals
   */


  double log_num, log_den, ln_acp, lnu, ld;
  /*
   * log_num: numerator for MH ratio
   * log_den: denumerator for MH ratio
   * ln_acp: MH ratio
   * lnu: random value for MH acceptance
   * ld: log determinant
   */

  arma::vec eta_p(dim);

  arma::mat loggamma_p(nobs, dim);
  loggamma_p = loggamma;
  //
  log_den = 0.0;
  it = 0;
  for(i = 0; i < nobs; i++){
    if(treatments(i) == t){
      if((curr_clu(it)) == (jj+1)){
        for(k = 0; k < dim; k++){
          log_den -= lgamma(exp(loggamma(i, k))) + exp(loggamma(i, k)) * log(JJ(i, k));
        }
      }
      it += 1;
    }
  }

  ld = logdet(sigma_star, dim);
  log_den += dmvnorm(eta, mu_star, sigma_star, dim, ld, 1);

  /*
   * propose new value for eta
   */

  for(k = 0; k < dim; k++){
    eta_p(k) = eta(k) + R::rnorm(0, mhtune);//R::runif(-1, 1);
  }

  it = 0;
  for(i = 0; i < nobs; i++){
    if(treatments(i) == t){
      if((curr_clu(it)-1) == (jj)){
        for(k = 0; k < dim; k++){
          loggamma_p(i, k) = loggamma(i, k) - eta(k) + eta_p(k);
        }
      }
      it += 1;
    }
  }

  log_num = 0.0;
  it = 0;
  for(i = 0; i < nobs; i++){
    if(treatments(i) == t){
      if((curr_clu(it)-1) == (jj)){
        for(k = 0; k < dim; k++){
          log_num -= lgamma(exp(loggamma_p(i, k))) + exp(loggamma_p(i, k)) * log(JJ(i, k));
        }
      }
      it += 1;
    }
  }
  log_num += dmvnorm(eta_p, mu_star, sigma_star, dim, ld, 1);

  ln_acp = log_num - log_den;

  lnu = log(R::runif(0.0, 1.0));

  if(lnu < ln_acp){
    // If accepted, update both eta and loggamma, and keep
    // track of acceptances
    eta_flag(jj) += 1;
    //Rcpp::Rcout << "accepted!" << std::endl;
    eta = eta_p;

    it = 0;
    for(i = 0; i < nobs; i++){
      if(treatments(i) == t){
        if((curr_clu(it)-1) == (jj)){
          loggamma.row(i) = loggamma_p.row(i);
        }
        it += 1;
      }
    }
  }//closes if accepted
  // Return output
  Rcpp::List eta_up(3);
  // eta, loggamma, acceptance,
  eta_up[0] = eta;
  eta_up[1] = loggamma;
  eta_up[2] = eta_flag;
  return eta_up;
}

/*
 * the following function updates the coefficients for prognostic covariates
 * $\boldsymbol{\beta}_{q}$
 */

Rcpp::List beta_update(arma::mat ZZ, arma::mat JJ, arma::mat loggamma,
                       arma::vec beta_temp, arma::mat beta_flag,
                       double mu_beta, arma::vec sigma_beta, int kk, double mhtune){

  // this function loops through K so it is called for one category at a time

  /*
   * ZZ: matrix of independent prognostic covariates
   * JJ: matrix of independent gamma variables
   * loggamma: matrix of log-linear predictor
   * beta_temp: vector of Q coefficient for kk-th category
   * beta_flag: flag for beta MH acceptance
   * mu_beta: prior mean
   * sigma_beta:  prior variance
   * kk: category we are working on
   * mhtune: tuning parameter for beta's proposal
   */
  int Q = ZZ.n_cols;
  int nobs = JJ.n_rows;
  int dim = JJ.n_cols;

  int i, q, h;
  /* indices for:
   * i: individuals
   * q: vovariates
   * h: instrumental to move along the vectorized matrix of beta coefficients
   */

  double log_num, log_den, ln_acp, lnu;
  /*
   * log_num: numerator for MH ratio
   * log_den: denumerator for MH ratio
   * ln_acp: MH ratio
   * lnu: random value for MH acceptance
   */

  double beta_p;

  arma::vec loggamma_p(nobs);

  for(q = 0; q < Q; q++){
    //h = q + kk * Q;
    h = kk + q * dim;

    log_den = 0.0;
    for(i = 0; i < nobs; i++){
        log_den = log_den - lgamma(exp(loggamma(i, kk))) + exp(loggamma(i, kk)) * log(JJ(i, kk));
    }
    log_den += R::dnorm4(beta_temp(q), mu_beta, sigma_beta(h), 1);

    // propose new value for beta (RW)
    beta_p = beta_temp(q) + R::rnorm(0, mhtune);

    for(i = 0; i < nobs; i++){
        loggamma_p(i) = loggamma(i, kk) - beta_temp(q) * ZZ(i, q) + beta_p * ZZ(i, q);
    }

    log_num = 0.0;
    for(i = 0; i < nobs; i++){
        log_num = log_num - lgamma(exp(loggamma_p(i))) + exp(loggamma_p(i)) * log(JJ(i, kk));
    }
    log_num += R::dnorm4(beta_p, mu_beta, sigma_beta(h), 1);

    ln_acp = log_num - log_den;
    /*if(ln_acp > 200.0){
      Rcpp::Rcout << "beta_p: " << beta_p << std::endl;
      Rcpp::Rcout << "loggamma: " << loggamma.col(kk).t() << std::endl;
      Rcpp::Rcout << "loggamma_p: " << loggamma_p.t() << std::endl;
      Rcpp::Rcout << "JJ: " << JJ.col(kk).t() << std::endl;
    }*/

    lnu = log(R::runif(0.0, 1.0));

    if(lnu < ln_acp){
      // If accepted, update both eta and loggamma, and keep track of acceptances
      beta_flag(q, kk) += 1;
      beta_temp(q) = beta_p;
      for(i = 0; i < nobs; i++){
        loggamma(i, kk) = loggamma_p(i);
      }
    }//closes if accepted
  }

  // Return output
  Rcpp::List beta_up(3);
  // eta, loggamma, acceptance,
  beta_up[0] = beta_temp;
  beta_up[1] = loggamma;
  beta_up[2] = beta_flag;
  return beta_up;
}


Rcpp::IntegerVector rmultinom_1(int size, Rcpp::NumericVector &probs, int N) {
  Rcpp::IntegerVector outcome(N);
  rmultinom(size, probs.begin(), N, outcome.begin());
  return outcome;
}

// [[Rcpp::export]]
arma::mat rmultinom_rcpp(int n, int size, arma::vec prob) {
  Rcpp::NumericVector probs(prob.size());
  probs = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(prob));
  unsigned int N = probs.length();
  arma::mat sim(N, n);
  Rcpp::NumericVector work(N);
  for (int i = 0; i < n; i++) {
    work = rmultinom_1(size, probs, N);
    sim.col(i) = Rcpp::as<arma::vec>(work);
  }
  return sim.t();
}

/*
 * Benchmark the function
 * adapted from https://gallery.rcpp.org/articles/recreating-rmultinom-and-rpois-with-rcpp/
 prob <- runif(200)
 prob <- prob/sum(prob) # standardize the probabilities
 size <- 500
 n <- 20

 set.seed(10)
 sim_r <- rmultinom(n, size, prob)
 set.seed(10)
 sim_rcpp <- rmultinom_rcpp(n, size, prob)
 all.equal(sim_r, sim_rcpp)

 microbenchmark::microbenchmark(
 rmultinom(1000, size, prob),
 rmultinom_rcpp(1000, size, prob)
 )
 */

// [[Rcpp::export]]
double dmultinom_rcpp(arma::vec x, int size, arma::vec prob, int Log){
  int K = prob.size();
  double s = arma::sum(prob);
  prob = prob/s;
  int N = 0;
  arma::vec i0(K, arma::fill::zeros);
  for(int i = 0; i < K; i++){
    x(i) = floor(x(i)+.5);
    N += x(i);
    if(prob(i) == 0){
      i0(i) = 1;
    }
  }

  /*if (any(i0)) {
   if (any(x[i0] != 0))
   return(if (log) -Inf else 0)
   if (all(i0))
   return(if (log) 0 else 1)
   x <- x[!i0]
   prob <- prob[!i0]
  }*/
  double r;
  r = lgamma(size + 1);
  for(int i = 0; i < K; i++){
    r += x(i) * log(prob(i)) - lgamma(x(i) + 1);
  }
  if (Log){
    return r;
  } else {
    return exp(r);
  }
}

arma::vec up_lambda_hs(arma::vec beta, arma::vec lambda, double tau){
  int len = beta.size();
  int k;
  double eta, upsi, tempps, ub, Fub, up;
  for(k = 0; k < len; k++){
    eta = pow(lambda(k), -2.0);
    upsi = R::runif(0, pow((1 + eta), - 1.0));
    tempps = pow(beta(k), 2.0)/(2*pow(tau, 2.0));
    ub = (1-upsi)/upsi;
    Fub = 1.0 - exp(-tempps * ub);
    if(Fub < pow(10.0, -4)){
      Fub = pow(10.0, -4);
    }
    up = R::runif(0, Fub);
    eta = -log(1.0 - up)/tempps;
    lambda(k) = pow(eta, -0.5);
  }
  return lambda;
}

double up_tau_hs(arma::vec beta, arma::vec lambda, double tau){
  int len = beta.size();
  double tempt, et, utau, ubt, Fubt, ut;
  tempt = arma::sum(pow((beta/lambda), 2.0))/2;
  et = pow(tau, -2.0);
  utau = R::runif(0.0, pow((1.0 + et), -1.0));
  ubt = (1.0 - utau)/utau;
  Fubt = R::pgamma(ubt, (((double) len) + 1.0)/2, pow(tempt, -1.0), 1, 0);
  Fubt = std::max(Fubt, pow(10.0, -4));
  ut = R::runif(0.0, Fubt);
  et = R::qgamma(ut, (((double) len) + 1.0)/2.0, pow(tempt, -1.0), 1, 0);
  if(et < .0001){
    et = .0001;
  }
  tau = pow(et, -0.5);
  if(tau < .0001){
    tau = .0001;
  }
  return tau;
}
