#pragma once

#include <Rcpp.h>

namespace Util {
  inline int rbernoulli(double p) {
    // NOTE: More efficient than R, but not equivalent.
    return (R::unif_rand() > p) ? 0 : 1;
    // Equivalent to R:
    //  double pp = std::max(p, 1-p);
    //  int ans = R::unif_rand() < pp ? 0 : 1;
    //  return (p > 0.5) ? 1-ans : ans;
  }
  inline double rgamma(double shape, double rate) { return R::rgamma(shape, 1/rate); }
  Rcpp::NumericVector rmvnorm(Rcpp::NumericVector mu, Rcpp::NumericMatrix eig_sigma);

  /* Compute the log likelihood ratio for a poisson with rates lambda1 and lambda0
   * lambda1 is the proposal, lambda0 is the current
   */
  inline double loglik_pois(int y, int n, double loglambda0, double loglambda1) {
    return y * (loglambda1 - loglambda0) - n * (::exp(loglambda1) - ::exp(loglambda0));
  }
}
