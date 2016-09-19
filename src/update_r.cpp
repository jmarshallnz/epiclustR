#include <Rcpp.h>
#include "util.h"

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
                int i,
                List state,
                List prior) {

  RNGScope scope;

  // extract information from our state
  NumericVector R = state["R"];

  double aR = prior["aR"];
  double bR = prior["bR"];

  // Gibbs update for kR
  double kR = Util::rgamma(aR + 0.5*(R.length()-2), bR+0.5*sum_r_squared(R));
  state["kR"] = kR;

  return state;
}
