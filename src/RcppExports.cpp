// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// myppmx
Rcpp::List myppmx(int iter, int burn, int thin, int nobs, int ncon, int ncat, arma::vec catvec, double alpha, int maug, int reuse, int cohesion, int similarity, int consim, arma::vec y, arma::vec xcon, arma::vec xcat, int npred, arma::mat xconp, arma::mat xcatp, arma::vec similparam, arma::vec modelpriors, arma::vec mhtune, int calibration);
RcppExport SEXP _treatppmx_myppmx(SEXP iterSEXP, SEXP burnSEXP, SEXP thinSEXP, SEXP nobsSEXP, SEXP nconSEXP, SEXP ncatSEXP, SEXP catvecSEXP, SEXP alphaSEXP, SEXP maugSEXP, SEXP reuseSEXP, SEXP cohesionSEXP, SEXP similaritySEXP, SEXP consimSEXP, SEXP ySEXP, SEXP xconSEXP, SEXP xcatSEXP, SEXP npredSEXP, SEXP xconpSEXP, SEXP xcatpSEXP, SEXP similparamSEXP, SEXP modelpriorsSEXP, SEXP mhtuneSEXP, SEXP calibrationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< int >::type ncon(nconSEXP);
    Rcpp::traits::input_parameter< int >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type catvec(catvecSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type maug(maugSEXP);
    Rcpp::traits::input_parameter< int >::type reuse(reuseSEXP);
    Rcpp::traits::input_parameter< int >::type cohesion(cohesionSEXP);
    Rcpp::traits::input_parameter< int >::type similarity(similaritySEXP);
    Rcpp::traits::input_parameter< int >::type consim(consimSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xcon(xconSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xcat(xcatSEXP);
    Rcpp::traits::input_parameter< int >::type npred(npredSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xconp(xconpSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xcatp(xcatpSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type similparam(similparamSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type modelpriors(modelpriorsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mhtune(mhtuneSEXP);
    Rcpp::traits::input_parameter< int >::type calibration(calibrationSEXP);
    rcpp_result_gen = Rcpp::wrap(myppmx(iter, burn, thin, nobs, ncon, ncat, catvec, alpha, maug, reuse, cohesion, similarity, consim, y, xcon, xcat, npred, xconp, xcatp, similparam, modelpriors, mhtune, calibration));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_treatppmx_myppmx", (DL_FUNC) &_treatppmx_myppmx, 23},
    {NULL, NULL, 0}
};

RcppExport void R_init_treatppmx(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
