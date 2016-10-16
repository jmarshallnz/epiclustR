#include "expected_cases.h"

using namespace Rcpp;

NumericMatrix log_case_rate(const Data &data, List state, bool smoothed, IntegerVector urange) {

  double       fe = state["fe"];
  NumericVector R = state["R"];
  NumericVector U = state["U"];
  IntegerMatrix X = state["X"];
  NumericVector betaX = state["betaX"];

  NumericMatrix e = no_init(data.cases.nrow(), urange.length());
  for (int t = 0; t < e.nrow(); t++) {
    for (int i = 0; i < urange.length(); i++) {
      int u = urange[i]-1;
      double log_e = fe + R[t] + U[u];
      if (!smoothed)
        log_e += X(t,data.mbrg[u]-1) * betaX[data.mbrg[u]-1];
      e(t, i) = log_e;
    }
  }
  return e;
}

//' Computed expected cases per time per space
//'
//' @param data a list containing the data
//' @param state the current state of the Markov chain
//' @param smoothed set true to smooth the cases (i.e. preclude outbreaks). Defaults to false
//' @param spatial the spatial locations to compute the expected cases at.
//' @return a matrix containing the expected cases in time and space
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix expected_cases(List data, List state, bool smoothed = false, IntegerVector spatial = IntegerVector(0)) {

  // setup our objects
  Data d(data);

  if (spatial.length() == 0) {
    spatial = seq_along(d.mbrg);
  }

  NumericMatrix log_lambda = log_case_rate(d, state, smoothed, spatial);

  NumericMatrix e = no_init(log_lambda.nrow(), log_lambda.ncol());
  for (int t = 0; t < e.nrow(); t++) {
    for (int u = 0; u < e.ncol(); u++) {
      e(t, u) = d.n[spatial[u]-1] * ::exp(log_lambda(t,u));
    }
  }

  return e;
}

//' Compute expected cases over time
//'
//' @param data a list containing the data
//' @param state the current state of the Markov chain
//' @param smoothed set true to smooth the cases (i.e. preclude outbreaks). Defaults to false
//' @param spatial the set of spatial regions to compute the expected cases at
//' @return a matrix containing the expected cases in time and space
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector cases_per_time(List data, List state, bool smoothed = false, IntegerVector spatial = IntegerVector(0)) {
  NumericMatrix cases = expected_cases(data, state, smoothed, spatial);
  NumericVector out = no_init(cases.nrow());
  for (int i = 0; i < cases.nrow(); i++) {
    out[i] = sum(cases(i,_));
  }
  return out;
}
