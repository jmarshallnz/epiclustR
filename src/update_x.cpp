#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

inline int rbernoulli(double p) {
  double pp = std::max(p, 1-p);
  int ans = R::unif_rand() < pp ? 0 : 1;
  return (p > 0.5) ? 1-ans : ans;

  // a more efficient implementation (but different to R's) is just. Saves a few compares
//  return (R::unif_rand() > p) ? 0 : 1;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix x_sample_rux2(NumericMatrix cases,
                                  NumericVector n,
                                  double fe,
                                  NumericVector R,
                                  NumericVector U,
                                  NumericVector betaX,
                                  double pX,
                                  List rgmb) {
  RNGScope scope;
  const int n_t = R.length();
  const int n_r = rgmb.length();
  IntegerMatrix X = no_init(n_t, n_r);
  for (int r = 0; r < n_r; r++) {
    for (int t = 0; t < n_t; t++) {
      double loglr = 0;
      NumericVector rgmb_r = rgmb[r];
      for (int j = 0; j < rgmb_r.length(); j++) {
        int u = rgmb_r[j]-1;
        double loglambda0 = fe + R[t] + U[u];
        double loglambda1 = loglambda0 + betaX[r];
        loglr += n[u] * (exp(loglambda1) - exp(loglambda0)) - cases(t,u) * betaX[r];
      }
      double p = pX / (exp(loglr) * (1-pX) + pX);
      X(t, r) = rbernoulli(p);
    }
  }
  return X;
//  lenR <- length(R)
//  lambda0 <- rep(n,each=lenR)*exp(fe+rep(R,ncol(n))+rep(U,each=lenR))
//  lambda1 <- lambda0 * exp(rep(betaX[mbrg],each=lenR))
//  squashProd(cases * (0 - rep(betaX[mbrg],each=lenR)) - lambda0 + lambda1)
}
