#include "util.h"

using namespace Rcpp;

namespace Util {
  Rcpp::NumericVector rmvnorm(NumericVector mu, NumericMatrix eig_sigma) {
    NumericVector Z = no_init(mu.length());
    for (int i = 0; i < Z.length(); i++)
      Z[i] = R::norm_rand();

    NumericVector X = mu;
    for (int i = 0; i < X.length(); i++) {
      for (int j = 0; j < X.length(); j++) {
        X[i] += eig_sigma(i,j) * Z[j];
      }
    }

    return X;
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector rmvnorm(NumericVector mu, NumericMatrix eig_sigma) {
// X <- rnorm(length(mu));
// as.numeric(mu + eigen_sigma %*% X)
  RNGScope scope;
  return Util::rmvnorm(mu, eig_sigma);
}
