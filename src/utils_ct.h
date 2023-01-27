double quform(arma::vec, arma::vec, int);
double inner_product(arma::vec, int, arma::vec, int, int);
double squared_norm(arma::vec, int, int, int);
arma::vec cholesky(arma::vec, int);
double logdet(arma::vec, int);
double dinvgamma(double, double, double, int);
double dN_IG(double, double, double, double, double, double, int);
double dmvnorm(arma::vec, arma::vec, arma::vec, int,  double,  int);
arma::vec ran_mvnorm(arma::vec, arma::vec, int);
arma::vec dmvnrm(arma::mat, arma::rowvec, arma::mat, bool);
double dinvwish(arma::vec, int, double, double, int, int);
arma::vec ran_iwish(int, arma::vec, int);
double gsimconNN(double, double, double, double, double, int, int, int);
double gsimconNNIG(double, double, double, double, double,
                   double, int, int, int);
double gsimcatDM(arma::vec, arma::vec, int, int, int);
//Rcpp::List ranppmx(int, int, int, double, int, int, arma::vec, arma::vec, arma::vec, double,
//      double, double, double, double, arma::vec);
double calculate_gamma(arma::mat, arma::mat, arma::vec, int, int, int, int);
//Rcpp::List eta_update(arma::mat, arma::mat, arma::vec, //arma::vec,
//                      arma::vec, int,
//                      arma::vec, arma::vec, arma::vec, arma::vec, int, double);
Rcpp::List eta_update(arma::mat, arma::mat, arma::vec, arma::vec, int,
                       arma::vec, arma::vec,
                       arma::vec, arma::mat, int, double);

Rcpp::List beta_update(arma::mat, arma::mat, arma::mat, arma::vec, arma::mat,
                       double, arma::vec, int, double);
double log_mult(arma::mat, arma::mat);
Rcpp::IntegerVector rmultinom_1(int, Rcpp::NumericVector, int);
arma::mat rmultinom_rcpp(int, int, arma::vec);
double dmultinom_rcpp(arma::vec, int, arma::vec, int);
double myround(double);
arma::vec up_lambda_hs(arma::vec, arma::vec, double);
double up_tau_hs(arma::vec, arma::vec, double);
