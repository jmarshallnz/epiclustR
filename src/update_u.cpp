#include <Rcpp.h>
#include "util.h"

using namespace Rcpp;

double u_likelihood(NumericMatrix cases,
                         NumericVector n,
                         double fe,
                         NumericVector R,
                         IntegerMatrix X,
                         NumericVector mbrg,
                         NumericVector betaX,
                         double curr,
                         double prop,
                         int u) {
  const int n_r = R.length();
  double lr = 0;
  for (int t = 0; t < n_r; t++) {
    double loglambda = fe + R[t] + betaX[mbrg[u]-1] * X(t, mbrg[u]-1);
    lr += cases(t,u) * (prop - curr)
          - n[u] * (exp(loglambda+prop) - exp(loglambda+curr));
  }
  return lr;
//  tps <- length(R)
//  lambda_curr <- n[j]*exp(fe+R+curr+rep(betaX[mbrg[j]],tps)*X[,mbrg[j]])
//  lambda_prop <- lambda_curr * exp(prop-curr)
//  sum(cases[,j] * (prop - curr) - lambda_prop + lambda_curr)
}

// [[Rcpp::export]]
double sum_u_squared(NumericVector U, NumericMatrix nb) {
  double sum = 0;
  for (int i = 0; i < U.length(); i++) {
    for (int j = 1; j < nb(i,0)+1; j++) {
      double d = U[i] - U[nb(i, j)-1];
      sum += d*d;
    }
  }
  return 0.5*sum; // We divide by 2 here as the above sum will count all squared distances twice.
}

// [[Rcpp::export]]
Rcpp::List update_u(NumericMatrix cases,
                NumericVector n,
                NumericVector mbrg,
                NumericMatrix nb,
                int i,
                List state,
                List prior) {

  RNGScope scope;

  // save our state
  List out = state;
  double       fe = state["fe"];
  NumericVector R = state["R"];
  NumericVector U = state["U"];
  IntegerMatrix X = state["X"];
  NumericVector betaX = state["betaX"];
  double kU = state["kU"];
  NumericVector acceptU = state["acceptU"];
  NumericVector rejectU = state["rejectU"];

  double aU = prior["aU"];
  double bU = prior["bU"];
  double sigmaU = prior["sigmaU"];

  // Gibb's step to update kU
  kU = Util::rgamma(aU + 0.5*(U.length()-1), bU + 0.5*sum_u_squared(U, nb));
  out["kU"] = kU;

  for (int u = 0; u < U.length(); u++) {
    double proposal, ap;
    if (i % 2 == 0) { // propose from the prior
      double sumU = 0;
      for (int j = 1; j < nb(u,0) + 1; j++)
        sumU += U[nb(u,j)-1];

      proposal = R::rnorm(sumU/nb(u,0), 1/sqrt(kU * nb(u,0)));
      ap = u_likelihood(cases, n, fe, R, X, mbrg, betaX, U[u], proposal, u);
    } else { // M-H random walk proposal
      proposal = R::rnorm(U[u], sigmaU);
      double sumU = 0;
      for (int j = 1; j < nb(u,0) + 1; j++) {
        double d1 = U[nb(u,j)-1] - proposal;
        double d2 = U[nb(u,j)-1] - U[u];
        sumU += d1*d1 - d2*d2;
      }
      ap = u_likelihood(cases, n, fe, R, X, mbrg, betaX, U[u], proposal, u) - kU * sumU / 2;
    }
    double un = R::unif_rand();
    if (ap >= 0 || un<=exp(ap)) {
      U[u] = proposal;
      acceptU[i%2]++;
    } else {
      rejectU[i%2]++;
    }
  }
  double meanU = mean(U);
  out["fe"] = fe + meanU;
  out["U"] = U - meanU;
  out["acceptU"] = acceptU;
  out["rejectU"] = rejectU;

  return out;
}
