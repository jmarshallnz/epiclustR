#include "update.h"
#include "util.h"

using namespace Rcpp;

double u_likelihood(const Data &data,
                    const State &s,
                         double prop,
                         int u) {
  const int n_r = s.R.length();
  double lr = 0;
  for (int t = 0; t < n_r; t++) {
    double loglambda = s.fe + s.R[t] + s.betaX[data.mbrg[u]-1] * s.X(t, data.mbrg[u]-1);
    lr += Util::loglik_pois(data.cases(t,u), data.n[u], loglambda+s.U[u], loglambda+prop);
  }
  return lr;
//  tps <- length(R)
//  lambda_curr <- n[j]*exp(fe+R+curr+rep(betaX[mbrg[j]],tps)*X[,mbrg[j]])
//  lambda_prop <- lambda_curr * exp(prop-curr)
//  sum(cases[,j] * (prop - curr) - lambda_prop + lambda_curr)
}

double sum_u_squared(const NumericVector &U, const NumericMatrix &nb) {
  double sum = 0;
  for (int i = 0; i < U.length(); i++) {
    for (int j = 1; j < nb(i,0)+1; j++) {
      double d = U[i] - U[nb(i, j)-1];
      sum += d*d;
    }
  }
  return 0.5*sum; // We divide by 2 here as the above sum will count all squared distances twice.
}

void update_u(const Data &data,
                int i,
                State &s,
                List prior,
                List control) {

  RNGScope scope;

  double aU = prior["aU"];
  double bU = prior["bU"];
  double sigmaU = control["sigmaU"];

  // Gibb's step to update kU
  s.kU = Util::rgamma(aU + 0.5*(s.U.length()-1), bU + 0.5*sum_u_squared(s.U, data.nb));

  for (int u = 0; u < s.U.length(); u++) {
    double proposal, ap;
    if (i % 2 == 0) { // propose from the prior
      double sumU = 0;
      for (int j = 1; j < data.nb(u,0) + 1; j++)
        sumU += s.U[data.nb(u,j)-1];

      proposal = R::rnorm(sumU/data.nb(u,0), 1/::sqrt(s.kU * data.nb(u,0)));
      ap = u_likelihood(data, s, proposal, u);
    } else { // M-H random walk proposal
      proposal = R::rnorm(s.U[u], sigmaU);
      double sumU = 0;
      for (int j = 1; j < data.nb(u,0) + 1; j++) {
        double d1 = s.U[data.nb(u,j)-1] - proposal;
        double d2 = s.U[data.nb(u,j)-1] - s.U[u];
        sumU += d1*d1 - d2*d2;
      }
      ap = u_likelihood(data, s, proposal, u) - s.kU * sumU / 2;
    }
    double un = R::unif_rand();
    if (ap >= 0 || un <= ::exp(ap)) {
      s.U[u] = proposal;
      s.acceptU[i%2]++;
    } else {
      s.rejectU[i%2]++;
    }
  }
}
