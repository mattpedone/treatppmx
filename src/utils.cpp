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

  #include <RcppArmadillo.h>
  #include <math.h>
  #include <Rcpp.h>
  #include <R.h>
  #include <Rmath.h>
  
//quadratic form
double quform(arma::vec x, arma::vec A, int dim){
  
  int i, j;
  double sm = 0.0;
  for (i = 1; i < dim; i++){
    for (j = 0; j < i; j++){
      sm += x(i)*x(j)*A(i * dim + j);
    }
  }
  sm *= 2;
  for(i = 0; i < dim; i++){
    sm += x(i)*x(i)*A(i * dim + i);
  }
  return(sm);
}

//inner product
/* Each vector is of length n, and is stored in memory in successive locations 
 * that are at a given distance apart. For an ordinary vector, a distance of 1 
 * is used, but other distances can come about from looking at columns of a 
 * matrix.
 */
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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

//[[Rcpp::export]]
arma::mat wa(arma::vec v, int dim){
  arma::vec wv = v;
  arma::mat work(dim, dim);
  for(int row = 0; row < dim; row++){
    for(int col = 0; col < dim; col++){
      work(row, col) = wv(row * dim + col);
    }
  }
  
  work = arma::chol(work, "lower");
  //arma::log_det(ldSig, sign, work);
  return work;
}

//log determinant of a posdef sym matrix

/*
 * If the symmetric positive definite matrix A is represented by its Cholesky 
 * decomposition A = LL' or A = U'U, then the determinant of this matrix can 
 * be calculated as the product of squares of the diagonal elements of L or U.
 * vector A of dim (n x n) x 1 stores the cholesky decomposition of the 
 * (n x n) matrix whose log determinant we are looking for
 * A n x n matrix whose determinant we want to compute 
 */
// [[Rcpp::export]]
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

//[[Rcpp::export]]
double wa_det(arma::vec v, int dim){
  //arma::vec wv = v;
  arma::mat work(dim, dim);
  double res, sign;
  for(int row = 0; row < dim; row++){
    for(int col = 0; col < dim; col++){
      work(row, col) = v(row * dim + col);
    }
  }
  //work = arma::inv(work);
  //work = arma::chol(work, "lower");
  //work = work.t();
  //Rcpp::Rcout << work << std::endl;
  arma::log_det(res, sign, work);
  //res =log(det(work));
  return res;
}

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
  
  ldens = R::dnorm(mu, mu0, sqrt(sig2/k0),logout) + dinvgamma(sig2, a0, b0, logout);
  //	Rprintf("ldens = %f\n", ldens);
  if(logout){ return ldens;
  }else{return exp(ldens);}
  
}

/* Multivariate normal density
 * y is the observation vector for which the density will be computed
 * mu is mean vector
 * Sig è la VARIANCE matrix
 * dim dimension of multivariate distribution
 * ld vuole log determinant variance matrix
 * logout logical if 1 returns log density
 */
// [[Rcpp::export]]
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

// [[Rcpp::export]]
arma::vec ran_mvnorm_old(arma::vec m, arma::vec Sig, int dim){
  
  int i,j;
  arma::mat cholV(dim, dim);
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      cholV(i, j) = Sig(i * dim + j);
    }
  }
  cholV = arma::chol(cholV, "lower");
  arma::vec z(dim);
  arma::vec out(dim);
  for(i = 0; i < dim; i++){
    z(i) = R::rnorm(0,1);
    out(i) = m(i);
    for(j = 0; j <= i; j++){
      out(i) += cholV(i, j)*z(j);
    }
  }
  return out;
}

/* Sampling Multivariate Normal
 * m vector of dimension dim storing the mean
 * Sig vector for VARIANCE MATRIX
 * dim dimension of multivarite normal distribution
 */
// [[Rcpp::export]]
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
  // [[Rcpp::export]]
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

// [[Rcpp::export]]
arma::vec ran_iwish_old(int nu, arma::vec Sig, int dim){
  
  //Rcpp::Rcout << "input" << Sig << std::endl;
  
  int i, j, k;
  arma::mat cholS(dim, dim, arma::fill::zeros);
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      cholS(i, j) = Sig(i * dim + j);
    }
  }
  //Rcpp::Rcout << "chols" << cholS << std::endl;
  
  cholS = arma::chol(cholS, "lower");
  
  /* for(i = 0; i < dim; i++){
  *  for(j = 0; j < dim; j++){
  *    Sig(i * dim + j) = cholS(i, j);
  *  }
  }*/
  
  //Rcpp::Rcout << "Sig" << Sig << std::endl;
  arma::vec x(dim);
  arma::vec zeros(dim);
  zeros.fill(0.0);
  arma::vec outv(dim * dim);
  arma::mat out(dim, dim);
  out.fill(0.0);
  
  for(i = 0; i < nu; i++){
    x = ran_mvnorm(zeros, Sig, dim);
    for(j = 0; j < dim; j++){
      for(k = 0; k <= j; k++){
        out(j, k) += x(j)*x(k);
      }
    }
  }
  
  /* fill the upper triangular part with lower triangular part */
    for(j=0; j<dim; j++){
      for(k=0; k<j; k++){
        out(k, j) = out(j, k);
      }
    }
    
    out = arma::inv(out);
  
    for(j = 0; j < dim; j++){
      for(k = 0; k < dim; k++){
        outv(j * dim + k) = out(j, k);
      }
    }
    
  return outv;
}

/* Random Draw from Wishart distribution
 * nu degrees of freedom
 * Sig matrix parameter NOT THE CHOLESKY DECOMPOSITION
 * dim dimension of the Wishart distribution
 */
// [[Rcpp::export]]
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
double gsimconNN(double m0, double v2, double s20, double sumx, double sumx2, double mle,
                 int n,  int DD, int cal, int logout){
  
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
  
  ld5 = R::dnorm(mle, m0, sqrt(s20),1);
  ld6 = R::dnorm(mle, mus, sqrt(s2s),1);
  
  out = ld1 + ld2 - ld3;
  if(DD==1) out = ld1 + ld3 - ld4;
  if(cal==1) out = ld5 - ld6;
  if(!logout) out = exp(out);
  return(out);
  
}

/* normal-normal-IG Similarity function with x following normal and m,v a normal-IG.
 * I didn't carry out integration explicitly over m and v. I simply used the fact that
 * marginal likelihood (similarity) is equal to likelihood x prior / posterior.  This
 * requires inserting a value for m and v which I use 0 and 1.
 * The double dipper is included as an argument.
 */
double gsimconNNIG(double m0, double k0, double nu0, double s20, double sumx, double sumx2,
                   double mnmle, double s2mle, int n, int DD, int cal, int logout){
  
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
  
  ld5 = dN_IG(mnmle, s2mle, m0, k0, a0, b0, 1);
  ld6 = dN_IG(mnmle, s2mle, m0s, k0s, a0s, b0s, 1);
  
  out = ld1 + ld2 - ld3;

  if(DD==1) out = ld1 + ld3 - ld4;
  if(cal==1) out = ld5 - ld6;
  if(!logout) out = exp(out);
  
  return(out);
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
  
  return(out);
}

/* Multivariate normal-normal-inverse-wishart Similarity function with x following normal
 * and (m, V) following a N-IW.   I evaluate the ratio between likelihood multiplied by prior
 * divided by posterior.  Arguments are the following
 *
 * m0 - is the prior mean vector of dimension dim
 * lam0 - is the prior scale parameter that is a scalar
 * nu0 - is the prior degrees of freedom and is a scalar
 * V0 - is the prior covariance matrix of dimension dim x dim
 *
 * sumxvec - is a vector of dimension 1 x dim containing totals of covariates
 * SumSq - dim x dim matrix whose value is \sum_{i=1}^{n} (x_i - xbar) (x_i - xbar)'
 * sumsq - scalar whose value is \sum_{i=1}^{n} x_i'x_i
 * dim - is a scalar to provides dimension
 * n - is a scalar indication number of observations
 * DD - is a logical determining if double dipper should be used or auxiliary
 * logout - a logical indicating on whether log of g(x) should provided.
 *
 * The double dipper is included as an argument.
 *
 * NOTE: SINCE ARGUMENTS OF LIKELIHOOD ARE ARBITRARY, I USE THE FOLLOWING
 * m_j = zero vector
 * V_j = Identity
 * THIS IS REFLECTED IN THE VALUE OF THE "LIKELIHOOD"
 */
double gsimconMVN_MVNIW(arma::vec m0, double lam0, double nu0, arma::vec Sig0, 
                        int dim, arma::vec sumxvec, arma::vec SumSq, double sumsq, 
                        int n,  int DD, int logout){
  
  int i, j, jj;
  double ld1, ld2, ld3, ld4;
  double ldSig, sign, ldet0, ldets, ldetss;
  double out;
  arma::vec scr1(dim * dim);
  arma::mat work(dim, dim);
  arma::vec scr2(dim * dim); //*dim ?
  arma::vec scr3(dim);
  arma::vec scr4(dim * dim);
  arma::vec scr5(dim * dim);
  double lams, lamss, nus, nuss;
  
  // Loglikelihood evaluation
  ld1 = -0.5*n*dim*log(2*M_PI) -  0.5*sumsq;
  
  // This is the log prior evaluation for auxiliary
  // Note, the first argument of dmvnorm is m_j.  Since this value can be selected
  // arbitrarily I set it to m0
  for(j = 0; j < dim * dim; j++){ 
      scr1(j) = (lam0)*Sig0(j);
  } // scr1 = ((1/lam0)Sig_j)^{-1}
  ldet0 = -dim*log(lam0);
  //qui va messa la vvarianza invece della precisione (scr1) e ldet0 va messo con segno +
  ld2 = dmvnorm(m0, m0, scr1, dim, ldet0, 1) + dinvwish(Sig0, dim, 1.0,  1.0, nu0, 1);
  
  // Auxiliary
  lams = lam0 + n;
  nus = nu0 + n;
  
  for(j = 0; j < dim; j++){
    scr1(j) = (1.0/n)*sumxvec(j) - m0(j); // scr1 is now (xbar-m0)
  }
  
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      work(i, j) = scr1(i * dim +j);
    }
  }
  //qui mi serve matrix product per fare tra 2 vettori il prodotto tra matrici
  work = work * work;//scr1 * scr1;//matrix_product(scr1, scr1, scr2, dim, dim, 1); // scr2 is (xbar-m0)(xbar-m0)'
  
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
       scr2(i * dim +j) = work(i, j);
    }
  }
  for(j = 0; j < dim; j++){
    scr3(j) = (lam0*m0(j) + sumxvec(j))/(lam0 + n); // scr3 is Ms
      for(jj = 0; jj < dim; jj++){
      scr4(j * dim + jj) = Sig0(j * dim + jj) + SumSq(j * dim + jj) + ((lam0*n)/(lam0 + n))*scr2(j *dim + jj); // scr4 is Psis
      scr5(j * dim + jj) = scr4(j * dim + jj);
      scr1(j * dim + jj) = (lams)*Sig0(j * dim + jj); // Now scr1 = ((1/lams)Sig_j)^{-1}
      }
  }
  ldets = -dim*log(lams);
  
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      work(i, j) = scr5(i * dim +j);
    }
  }
  
  work = arma::chol(work, "lower");
  
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      scr5(i * dim +j) = work(i, j);
    }
  }
  arma::log_det(ldSig, sign, work);
  
  // Note, the first argument is m_j in notes which is arbitrary so I set it to m0
  //vd commento prima x dmvnorm
  ld3 = dmvnorm(m0, scr3, scr1, dim, ldets, 1) + dinvwish(scr4, dim, 1.0, exp(ldSig*sign), nus, 1);
  
  // double dipper
  lamss = lams + n;
  nuss = nus + n;
  
  for(j = 0; j < dim; j++){
    scr1(j) = (1.0/n)*sumxvec(j) - scr3(j); // scr1 = (xbar-Ms) recall scr3=Ms
  }
  
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      work(i, j) = scr1(i * dim +j);
    }
  }
  //qui mi serve matrix product per fare tra 2 vettori il prodotto tra matrici
  work = work * work;//scr1 * scr1;//matrix_product(scr1, scr1, scr2, dim, dim, 1); // scr2 is (xbar-m0)(xbar-m0)'
  
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      scr2(i * dim +j) = work(i, j);
    }
  }
  
  //scr2 = scr1 * scr1;//matrix_product(scr1, scr1, scr2, dim, dim, 1); // scr2 = (xbar-Ms)(xbar-Ms)'
  
  for(j = 0; j < dim; j++){
    scr1(j) = (lams*scr3(j) + sumxvec(j))/(lams + n); // scr1 = Mss
    for(jj = 0; jj < dim; jj++){
      
      scr5(j * dim + jj) = scr4(j * dim + jj) + SumSq(j * dim + jj) +  // scr4 = Psis, scr5=Psiss
        ((lams*n)/(lams + n))*scr2(j * dim + jj); //
      
    }
  }
  
  //CONTROLLA!!!! scr 3 ora mi serve come matrice (posso usare scr4)
  for(j = 0; j < dim * dim; j++){
    scr3(j) = scr5(j); // extra Psiss to get log determinant
    scr2(j) = (lamss)*Sig0(j); // scr2 = ((1/lamss)Sig_j)^{-1}
    }
  ldetss = -dim*log(lamss);

    for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      work(i, j) = scr4(i * dim +j);
    }
  }
  
  work = arma::chol(work, "lower");
  arma::log_det(ldSig, sign, work);
  
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      scr4(i * dim +j) = work(i, j);
    }
  }
  
  // Note, the first argument is m_j in notes which is arbitrary so I set it to m0
  //vd commento prima x dmvnorm
  ld4 = dmvnorm(m0, scr1, scr2, dim, ldetss, 1) + dinvwish(scr5, dim, 1.0, exp(ldSig*sign), nuss, 1);
  
  out = ld1 + ld2 - ld3;
  
  if(DD==1) out = ld1 + ld3 - ld4;
  if(!logout) out = exp(out);
  //	Rprintf("out = %f\n", out);
  return(out);
  
}
