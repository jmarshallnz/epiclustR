#include "util.h"

using namespace Rcpp;

// [[Rcpp::export]]
int rbernoulli(double p) {
  RNGScope scope;
  return Util::rbernoulli(p);
}
