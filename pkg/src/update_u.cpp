#include "update.h"
#include "util.h"

using namespace Rcpp;

double u_likelihood(const Data &data,
                         double fe,
                         NumericVector R,
                         IntegerMatrix X,
                         NumericVector betaX,
                         double curr,
                         double prop,
                         int u) {
  const int n_r = R.length();
  double lr = 0;
  for (int t = 0; t < n_r; t++) {
    double loglambda = fe + R[t] + betaX[data.mbrg[u]-1] * X(t, data.mbrg[u]-1);
    lr += data.cases(t,u) * (prop - curr)
          - data.n[u] * (::exp(loglambda+prop) - ::exp(loglambda+curr));
  }
  return lr;
//  tps <- length(R)
//  lambda_curr <- n[j]*exp(fe+R+curr+rep(betaX[mbrg[j]],tps)*X[,mbrg[j]])
//  lambda_prop <- lambda_curr * exp(prop-curr)
//  sum(cases[,j] * (prop - curr) - lambda_prop + lambda_curr)
}

double sum_u_squared(NumericVector U, const NumericMatrix &nb) {
  double sum = 0;
  for (int i = 0; i < U.length(); i++) {
    for (int j = 1; j < nb(i,0)+1; j++) {
      double d = U[i] - U[nb(i, j)-1];
      sum += d*d;
    }
  }
  return 0.5*sum; // We divide by 2 here as the above sum will count all squared distances twice.
}

Rcpp::List update_u(const Data &data,
                int i,
                List state,
                List prior,
                List control) {

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
  double sigmaU = control["sigmaU"];

  // Gibb's step to update kU
  kU = Util::rgamma(aU + 0.5*(U.length()-1), bU + 0.5*sum_u_squared(U, data.nb));
  out["kU"] = kU;

  for (int u = 0; u < U.length(); u++) {
    double proposal, ap;
    if (i % 2 == 0) { // propose from the prior
      double sumU = 0;
      for (int j = 1; j < data.nb(u,0) + 1; j++)
        sumU += U[data.nb(u,j)-1];

      proposal = R::rnorm(sumU/data.nb(u,0), 1/::sqrt(kU * data.nb(u,0)));
      ap = u_likelihood(data, fe, R, X, betaX, U[u], proposal, u);
    } else { // M-H random walk proposal
      proposal = R::rnorm(U[u], sigmaU);
      double sumU = 0;
      for (int j = 1; j < data.nb(u,0) + 1; j++) {
        double d1 = U[data.nb(u,j)-1] - proposal;
        double d2 = U[data.nb(u,j)-1] - U[u];
        sumU += d1*d1 - d2*d2;
      }
      ap = u_likelihood(data, fe, R, X, betaX, U[u], proposal, u) - kU * sumU / 2;
    }
    double un = R::unif_rand();
    if (ap >= 0 || un <= ::exp(ap)) {
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
