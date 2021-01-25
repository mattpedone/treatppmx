double dinvgamma(double, double, double, int);
double dN_IG(double, double, double, double, double, double, int);
double dmvnorm(arma::vec, arma::vec, arma::vec, int,  double,  int); 
double gsimconNN(double, double, double, double, double, double, int, int, int, int);
double gsimconNNIG(double, double, double, double, double, double, double, 
                   double, int, int, int, int);
double gsimcatDM(arma::vec, arma::vec, int, int, int);
arma::vec ran_mvnorm(arma::vec, arma::vec, int); 
double dinvwish(arma::vec, int, double, double, int, int); 
arma::vec ran_iwish(int, arma::vec, int); 
double gsimconMVN_MVNIW(arma::vec, double, double, arma::vec, int, arma::vec, 
                        arma::vec, double, int, int, int); 