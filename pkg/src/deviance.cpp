#include "expected_cases.h"

using namespace Rcpp;

//' Deviance
//'
//' @param data a list containing the data
//' @param state the current state of the Markov chain
//' @return a matrix containing the log case rate in time and space
//' @export
// [[Rcpp::export]]
double deviance(List data, List state) {
  Data d(data);
  NumericMatrix log_lambda = log_case_rate(d, state);
  double dev = 0;
  for (int t = 0; t < log_lambda.nrow(); t++) {
    for (int u = 0; u < log_lambda.ncol(); u++) {
      dev += d.cases(t,u) * log_lambda(t,u) - d.n[u] * ::exp(log_lambda(t,u));
    }
  }
  return -2 * dev;
}
