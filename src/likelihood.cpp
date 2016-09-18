#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
double r_likelihood_rux2(NumericMatrix cases,
                         NumericVector n,
                         double fe,
                         NumericVector U,
                         NumericMatrix X,
                         NumericVector mbrg,
                         NumericVector betaX,
                         NumericVector curr,
                         NumericVector prop,
                         int j,
                         int k) {
  const int n_u = U.length();
  double lr = 0;
  for (int u = 0; u < n_u; u++) {
    double l_u = fe + U[u];
    for (int i = 0; i < curr.length(); i++) {
      double loglambda = l_u + X(i+j-1, mbrg[u]-1)*betaX[mbrg[u]-1];
      lr += cases(i+j-1,u) * (prop[i] - curr[i])
             - n[u] * (exp(loglambda + prop[i]) - exp(loglambda + curr[i]));
    }
  }
  return lr;
//  mbs <- ncol(n)
//   lambda_curr <- rep(n,each=k-j+1)*exp(fe+rep(curr,mbs)+rep(U,each=k-j+1)+X[j:k,mbrg]*rep(betaX[mbrg],each=k-j+1))
// lambda_prop <- lambda_curr * rep(exp(prop-curr),mbs)
// sum(cases[j:k,] * rep(prop - curr,mbs) - lambda_prop + lambda_curr)
}
