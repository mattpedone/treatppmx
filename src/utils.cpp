//define ARMA_DONT_PRINT_ERRORS 
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

//Multivariate normal density
// [[Rcpp::export]]
double dmvnorm(arma::vec y, arma::vec mu, arma::vec Sig, int dim, double ld, int logout){
  
  arma::vec scr(dim);
  int i;
  double qf, out;
  for(i = 0; i < dim; i++){
    scr(i) = y(i) - mu(i);
  }
  qf = quform(scr, Sig, dim);
  out = -(double) dim*M_LN_SQRT_2PI - 0.5*(ld + qf);
  if (logout)	return out;
  return exp(out);
}

//Sampling Multivariate Normal
// [[Rcpp::export]]
arma::vec ran_mvnorm(arma::vec m, arma::vec Sig, int dim){
  
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

/* The following provides a density function for the Inverse Wishart
distribution.  This function relies heavily on matrix.c
To use it I must create an matrix.o(object file)
Sig is argument and S is parameter.
*Sig - is SSigma^{-1} that is found in the trace of the density.  
*In Hoff is S_0\Sigma^{-1}
nu0 - degrees of freedom of the inverse-wishart function*/
  
  // [[Rcpp::export]]
double dinvwish(arma::vec Sig, int dim, double detSig, double detS, int nu0, 
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
      trace += Sig(i);
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

//Random Draw from Wishart distribution
// [[Rcpp::export]]
arma::mat ran_wish(int nu, arma::vec Sig, int dim){
  
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
  return out;
}

// Multivariate normal-normal-inverse-wishart Similarity function with x following normal
// and (m, V) following a N-IW.   I evaluate the ratio between likelihood multiplied by prior
// divided by posterior.  Arguments are the following
//
  // m0 - is the prior mean vector of dimension dim
// lam0 - is the prior scale parameter that is a scalar
// nu0 - is the prior degrees of freedom and is a scalar
// V0 - is the prior covariance matrix of dimension dim x dim
//
  // sumxvec - is a vector of dimension 1 x dim containing totals of covariates
// SumSq - dim x dim matrix whose value is \sum_{i=1}^{n} (x_i - xbar) (x_i - xbar)'
// sumsq - scalar whose value is \sum_{i=1}^{n} x_i'x_i
// dim - is a scalar to provides dimension
// n - is a scalar indication number of observations
// DD - is a logical determining if double dipper should be used or auxiliary
// logout - a logical indicating on whether log of g(x) should provided.
//
  // The double dipper is included as an argument.
//
  // NOTE: SINCE ARGUMENTS OF LIKELIHOOD ARE ARBITRARY, I USE THE FOLLOWING
// m_j = zero vector
// V_j = Identity
// THIS IS REFLECTED IN THE VALUE OF THE "LIKELIHOOD"

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
  //double *scr1, double *scr2, double *scr3, double *scr4, double *scr5
  
  double lams, lamss, nus, nuss;
  
  //	RprintVecAsMat("SumSq", SumSq, dim, dim);
  
  // Loglikelihood evaluation
  ld1 = -0.5*n*dim*log(2*M_PI) -  0.5*sumsq;
  
  // This is the log prior evaluation for auxiliary
  // Note, the first argument of dmvnorm is m_j.  Since this value can be selected
  // arbitrarily I set it to m0
  for(j = 0; j < dim * dim; j++){ 
      scr1(j) = (lam0)*Sig0(j);
  } // scr1 = ((1/lam0)Sig_j)^{-1}
  ldet0 = -dim*log(lam0);
  
  ld2 = dmvnorm(m0, m0, scr1, dim, ldet0, 1) + 
    dinvwish(Sig0, dim, 1.0,  1.0, nu0, 1);
  
  // Auxiliary
  lams = lam0 + n;
  nus = nu0 + n;
  
  for(j = 0; j < dim; j++){
    scr1(j) = (1.0/n)*sumxvec(j) - m0(j); // scr1 is now (xbar-m0)
  }
  
  
  //	RprintVecAsMat("xbar_m0", xbar_m0, 1, dim);
  
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      work(i, j) = scr1(i * dim +j);
    }
  }
  //qui mi serve matrix_product per fare tra 2 vettori il prodotto tra matrici
  work = work * work;//scr1 * scr1;//matrix_product(scr1, scr1, scr2, dim, dim, 1); // scr2 is (xbar-m0)(xbar-m0)'
  
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
       scr2(i * dim +j) = work(i, j);
    }
  }
  //	RprintVecAsMat("(xbar_m0) (xbar_m0)'", scr1, dim, dim);
  
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
  // scr1 = ((1/lams)Sig_j)^{-1}
  // scr3 = Ms
  // scr4 = Psis
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
  //qui mi serve matrix_product per fare tra 2 vettori il prodotto tra matrici
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
  
  // scr1 = Mss
  // scr2 = ((1/lamss)Sig_j)^{-1}
  // scr3 = not needed
  // scr4 = Psiss (used to get determinant
  // scr5 = Psiss
  
  //	Rprintf("dmvnorm(m0, Mss, lamssVj, dim, ldetss, scr2, 1) = %f\n", dmvnorm(m0, scr1, scr2, dim, ldetss, scr4, 1));
  
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
  //cholesky(scr3, dim, &ldV);
  
  // Note, the first argument is m_j in notes which is arbitrary so I set it to m0
  ld4 = dmvnorm(m0, scr1, scr2, dim, ldetss, 1) + dinvwish(scr5, dim, 1.0, exp(ldSig*sign), nuss, 1);
  
  out = ld1 + ld2 - ld3;
  
  if(DD==1) out = ld1 + ld3 - ld4;
  if(!logout) out = exp(out);
  //	Rprintf("out = %f\n", out);
  return(out);
  
}
