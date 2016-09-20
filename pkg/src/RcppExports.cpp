// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// expected_cases
Rcpp::NumericMatrix expected_cases(List data, List state, bool smoothed);
RcppExport SEXP epiclustR_expected_cases(SEXP dataSEXP, SEXP stateSEXP, SEXP smoothedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< List >::type state(stateSEXP);
    Rcpp::traits::input_parameter< bool >::type smoothed(smoothedSEXP);
    __result = Rcpp::wrap(expected_cases(data, state, smoothed));
    return __result;
END_RCPP
}
// cases_per_time
Rcpp::NumericVector cases_per_time(List data, List state, bool smoothed);
RcppExport SEXP epiclustR_cases_per_time(SEXP dataSEXP, SEXP stateSEXP, SEXP smoothedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< List >::type state(stateSEXP);
    Rcpp::traits::input_parameter< bool >::type smoothed(smoothedSEXP);
    __result = Rcpp::wrap(cases_per_time(data, state, smoothed));
    return __result;
END_RCPP
}
// rbernoulli
int rbernoulli(double p);
RcppExport SEXP epiclustR_rbernoulli(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    __result = Rcpp::wrap(rbernoulli(p));
    return __result;
END_RCPP
}
// rmvnorm
Rcpp::NumericVector rmvnorm(NumericVector mu, NumericMatrix eig_sigma);
RcppExport SEXP epiclustR_rmvnorm(SEXP muSEXP, SEXP eig_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eig_sigma(eig_sigmaSEXP);
    __result = Rcpp::wrap(rmvnorm(mu, eig_sigma));
    return __result;
END_RCPP
}
// update
Rcpp::List update(List data, int i, List state, List prior, List control);
RcppExport SEXP epiclustR_update(SEXP dataSEXP, SEXP iSEXP, SEXP stateSEXP, SEXP priorSEXP, SEXP controlSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< List >::type state(stateSEXP);
    Rcpp::traits::input_parameter< List >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< List >::type control(controlSEXP);
    __result = Rcpp::wrap(update(data, i, state, prior, control));
    return __result;
END_RCPP
}
