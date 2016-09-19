#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
double sum_r_squared(NumericVector R) {
  double sum = 0;
  for (int i = 0; i < R.length()-2; i++) {
    double r = R[i] - 2*R[i+1] + R[i+2];
    sum += r * r;
  }
  return sum;
}

// [[Rcpp::export]]
Rcpp::List update_r(NumericMatrix cases,
                NumericVector n,
                NumericVector mbrg,
                NumericMatrix nb,
                int i,
                List state,
                List prior) {

  RNGScope scope;

  return state;
}
