// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dm_ppmx
Rcpp::List dm_ppmx(int iter, int burn, int thin, int nobs, int PPMx, int ncon, int ncat, arma::vec catvec, double alpha, int CC, int reuse, int consim, int similarity, int calibration, arma::mat y, arma::vec xcon, arma::vec xcat, arma::vec similparam, arma::vec hP0_m0, arma::vec hP0_L0, double hP0_nu0, arma::vec hP0_V0, arma::vec mhtune);
RcppExport SEXP _treatppmx_dm_ppmx(SEXP iterSEXP, SEXP burnSEXP, SEXP thinSEXP, SEXP nobsSEXP, SEXP PPMxSEXP, SEXP nconSEXP, SEXP ncatSEXP, SEXP catvecSEXP, SEXP alphaSEXP, SEXP CCSEXP, SEXP reuseSEXP, SEXP consimSEXP, SEXP similaritySEXP, SEXP calibrationSEXP, SEXP ySEXP, SEXP xconSEXP, SEXP xcatSEXP, SEXP similparamSEXP, SEXP hP0_m0SEXP, SEXP hP0_L0SEXP, SEXP hP0_nu0SEXP, SEXP hP0_V0SEXP, SEXP mhtuneSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< int >::type PPMx(PPMxSEXP);
    Rcpp::traits::input_parameter< int >::type ncon(nconSEXP);
    Rcpp::traits::input_parameter< int >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type catvec(catvecSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type CC(CCSEXP);
    Rcpp::traits::input_parameter< int >::type reuse(reuseSEXP);
    Rcpp::traits::input_parameter< int >::type consim(consimSEXP);
    Rcpp::traits::input_parameter< int >::type similarity(similaritySEXP);
    Rcpp::traits::input_parameter< int >::type calibration(calibrationSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xcon(xconSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xcat(xcatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type similparam(similparamSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hP0_m0(hP0_m0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hP0_L0(hP0_L0SEXP);
    Rcpp::traits::input_parameter< double >::type hP0_nu0(hP0_nu0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hP0_V0(hP0_V0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mhtune(mhtuneSEXP);
    rcpp_result_gen = Rcpp::wrap(dm_ppmx(iter, burn, thin, nobs, PPMx, ncon, ncat, catvec, alpha, CC, reuse, consim, similarity, calibration, y, xcon, xcat, similparam, hP0_m0, hP0_L0, hP0_nu0, hP0_V0, mhtune));
    return rcpp_result_gen;
END_RCPP
}
// mvn_ppmx
Rcpp::List mvn_ppmx(int iter, int burn, int thin, int nobs, int PPMx, int ncon, int ncat, arma::vec catvec, double alpha, int CC, int reuse, int consim, int similarity, int calibration, arma::mat y, arma::vec xcon, arma::vec xcat, arma::vec similparam, arma::vec hP0_m0, arma::vec hP0_L0, double hP0_nu0, arma::vec hP0_V0, arma::vec mhtune);
RcppExport SEXP _treatppmx_mvn_ppmx(SEXP iterSEXP, SEXP burnSEXP, SEXP thinSEXP, SEXP nobsSEXP, SEXP PPMxSEXP, SEXP nconSEXP, SEXP ncatSEXP, SEXP catvecSEXP, SEXP alphaSEXP, SEXP CCSEXP, SEXP reuseSEXP, SEXP consimSEXP, SEXP similaritySEXP, SEXP calibrationSEXP, SEXP ySEXP, SEXP xconSEXP, SEXP xcatSEXP, SEXP similparamSEXP, SEXP hP0_m0SEXP, SEXP hP0_L0SEXP, SEXP hP0_nu0SEXP, SEXP hP0_V0SEXP, SEXP mhtuneSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< int >::type PPMx(PPMxSEXP);
    Rcpp::traits::input_parameter< int >::type ncon(nconSEXP);
    Rcpp::traits::input_parameter< int >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type catvec(catvecSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type CC(CCSEXP);
    Rcpp::traits::input_parameter< int >::type reuse(reuseSEXP);
    Rcpp::traits::input_parameter< int >::type consim(consimSEXP);
    Rcpp::traits::input_parameter< int >::type similarity(similaritySEXP);
    Rcpp::traits::input_parameter< int >::type calibration(calibrationSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xcon(xconSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xcat(xcatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type similparam(similparamSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hP0_m0(hP0_m0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hP0_L0(hP0_L0SEXP);
    Rcpp::traits::input_parameter< double >::type hP0_nu0(hP0_nu0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type hP0_V0(hP0_V0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mhtune(mhtuneSEXP);
    rcpp_result_gen = Rcpp::wrap(mvn_ppmx(iter, burn, thin, nobs, PPMx, ncon, ncat, catvec, alpha, CC, reuse, consim, similarity, calibration, y, xcon, xcat, similparam, hP0_m0, hP0_L0, hP0_nu0, hP0_V0, mhtune));
    return rcpp_result_gen;
END_RCPP
}
// myppmx
Rcpp::List myppmx(int iter, int burn, int thin, int nobs, int ncon, int ncat, arma::vec catvec, double alpha, int CC, int cohesion, int similarity, int consim, arma::vec y, arma::vec xcon, arma::vec xcat, arma::vec similparam, arma::vec modelpriors, arma::vec mhtune, int calibration);
RcppExport SEXP _treatppmx_myppmx(SEXP iterSEXP, SEXP burnSEXP, SEXP thinSEXP, SEXP nobsSEXP, SEXP nconSEXP, SEXP ncatSEXP, SEXP catvecSEXP, SEXP alphaSEXP, SEXP CCSEXP, SEXP cohesionSEXP, SEXP similaritySEXP, SEXP consimSEXP, SEXP ySEXP, SEXP xconSEXP, SEXP xcatSEXP, SEXP similparamSEXP, SEXP modelpriorsSEXP, SEXP mhtuneSEXP, SEXP calibrationSEXP) {
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
    Rcpp::traits::input_parameter< int >::type CC(CCSEXP);
    Rcpp::traits::input_parameter< int >::type cohesion(cohesionSEXP);
    Rcpp::traits::input_parameter< int >::type similarity(similaritySEXP);
    Rcpp::traits::input_parameter< int >::type consim(consimSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xcon(xconSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xcat(xcatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type similparam(similparamSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type modelpriors(modelpriorsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mhtune(mhtuneSEXP);
    Rcpp::traits::input_parameter< int >::type calibration(calibrationSEXP);
    rcpp_result_gen = Rcpp::wrap(myppmx(iter, burn, thin, nobs, ncon, ncat, catvec, alpha, CC, cohesion, similarity, consim, y, xcon, xcat, similparam, modelpriors, mhtune, calibration));
    return rcpp_result_gen;
END_RCPP
}
// calculate_gamma
double calculate_gamma(arma::mat eta, arma::vec curr_clu, int k, int i, int Log);
RcppExport SEXP _treatppmx_calculate_gamma(SEXP etaSEXP, SEXP curr_cluSEXP, SEXP kSEXP, SEXP iSEXP, SEXP LogSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type curr_clu(curr_cluSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type Log(LogSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_gamma(eta, curr_clu, k, i, Log));
    return rcpp_result_gen;
END_RCPP
}
// ranppmx
Rcpp::List ranppmx(int nobs, int similarity, int similparam, double alpha, int ncon, int ncat, arma::vec xcon, arma::vec xcat, arma::vec Cvec, double m0, double k0, double v0, double s20, double v, arma::vec dirweights);
RcppExport SEXP _treatppmx_ranppmx(SEXP nobsSEXP, SEXP similaritySEXP, SEXP similparamSEXP, SEXP alphaSEXP, SEXP nconSEXP, SEXP ncatSEXP, SEXP xconSEXP, SEXP xcatSEXP, SEXP CvecSEXP, SEXP m0SEXP, SEXP k0SEXP, SEXP v0SEXP, SEXP s20SEXP, SEXP vSEXP, SEXP dirweightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< int >::type similarity(similaritySEXP);
    Rcpp::traits::input_parameter< int >::type similparam(similparamSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type ncon(nconSEXP);
    Rcpp::traits::input_parameter< int >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xcon(xconSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xcat(xcatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Cvec(CvecSEXP);
    Rcpp::traits::input_parameter< double >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< double >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< double >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< double >::type s20(s20SEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dirweights(dirweightsSEXP);
    rcpp_result_gen = Rcpp::wrap(ranppmx(nobs, similarity, similparam, alpha, ncon, ncat, xcon, xcat, Cvec, m0, k0, v0, s20, v, dirweights));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_treatppmx_dm_ppmx", (DL_FUNC) &_treatppmx_dm_ppmx, 23},
    {"_treatppmx_mvn_ppmx", (DL_FUNC) &_treatppmx_mvn_ppmx, 23},
    {"_treatppmx_myppmx", (DL_FUNC) &_treatppmx_myppmx, 19},
    {"_treatppmx_calculate_gamma", (DL_FUNC) &_treatppmx_calculate_gamma, 5},
    {"_treatppmx_ranppmx", (DL_FUNC) &_treatppmx_ranppmx, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_treatppmx(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
