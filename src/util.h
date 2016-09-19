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
}
