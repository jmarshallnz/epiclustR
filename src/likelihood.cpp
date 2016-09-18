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

// [[Rcpp::export]]
double u_likelihood_rux2(NumericMatrix cases,
                         NumericVector n,
                         double fe,
                         NumericVector R,
                         NumericMatrix X,
                         NumericVector mbrg,
                         NumericVector betaX,
                         double curr,
                         double prop,
                         int j) {
  const int n_r = R.length();
  double lr = 0;
  for (int t = 0; t < n_r; t++) {
    double loglambda = fe + R[t] + betaX[mbrg[j-1]-1] * X(t, mbrg[j-1]-1);
    lr += cases(t,j-1) * (prop - curr)
          - n[j-1] * (exp(loglambda+prop) - exp(loglambda+curr));
  }
  return lr;
//  tps <- length(R)
//  lambda_curr <- n[j]*exp(fe+R+curr+rep(betaX[mbrg[j]],tps)*X[,mbrg[j]])
//  lambda_prop <- lambda_curr * exp(prop-curr)
//  sum(cases[,j] * (prop - curr) - lambda_prop + lambda_curr)
}

// [[Rcpp::export]]
double betax_likelihood_rux2(NumericMatrix cases,
                             NumericVector n,
                             double fe,
                             NumericVector R,
                             NumericVector U,
                             NumericMatrix X,
                             NumericVector rgmb_j,
                             double curr,
                             double prop,
                             int j) {
  const int n_t = R.length();
  const int n_r = rgmb_j.length();
  double lr = 0;
  for (int t = 0; t < n_t; t++) {
    for (int r = 0; r < n_r; r++) {
      int u = rgmb_j[r]-1;
      double loglambda = fe + R[t] + U[u];
      lr += cases(t,u) * X(t,j-1) * (prop - curr)
        - n[u] * (exp(loglambda + X(t,j-1)*prop) - exp(loglambda + X(t,j-1)*curr));
    }
  }
  return lr;
  // Xj <- X[rep(j-1,each=tps)*tps+rep(1:tps,lwch[j])]
  // lambda_curr <- rep(n[wch[[j]]],each=tps)*exp(fe+rep(R,lwch[j])+rep(U[wch[[j]]],each=tps)+Xj*curr)
  // lambda_prop <- lambda_curr * exp(Xj*(prop-curr))
  // sum(cases[,wch[[j]]] * X * (prop - curr) - lambda_prop + lambda_curr)
}
