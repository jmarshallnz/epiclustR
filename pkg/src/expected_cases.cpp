#include "expected_cases.h"

using namespace Rcpp;

NumericMatrix log_case_rate(const Data &data, List state, bool smoothed) {

  double       fe = state["fe"];
  NumericVector R = state["R"];
  NumericVector U = state["U"];
  IntegerMatrix X = state["X"];
  NumericVector betaX = state["betaX"];

  NumericMatrix e = no_init(data.cases.nrow(), data.cases.ncol());
  for (int t = 0; t < e.nrow(); t++) {
    for (int u = 0; u < e.ncol(); u++) {
      double log_e = fe + R[t] + U[u];
      if (!smoothed)
        log_e += X(t,data.mbrg[u]-1) * betaX[data.mbrg[u]-1];
      e(t, u) = log_e;
    }
  }
  return e;
}

//' Computed expected cases per time per space
//'
//' @param data a list containing the data
//' @param state the current state of the Markov chain
//' @param smoothed set true to smooth the cases (i.e. preclude outbreaks). Defaults to false
//' @return a matrix containing the expected cases in time and space
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix expected_cases(List data, List state, bool smoothed = false) {

  // setup our objects
  Data d(data);
  NumericMatrix log_lambda = log_case_rate(d, state, smoothed);

  NumericMatrix e = no_init(log_lambda.nrow(), log_lambda.ncol());
  for (int t = 0; t < e.nrow(); t++) {
    for (int u = 0; u < e.ncol(); u++) {
      e(t, u) = d.n[u] * ::exp(log_lambda(t,u));
    }
  }

  return e;
}
