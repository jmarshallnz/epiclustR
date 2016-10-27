#include "update.h"
#include "util.h"

using namespace Rcpp;

double u_likelihood(const Data &data,
                    const State &s,
                         double prop,
                         int u, int p) {
  double lr = 0;
  const NumericVector &period_p = data.p2t[p];
  const int n_p = period_p.length();
  for (int j = 0; j < n_p; j++) {
    const int t = period_p[j]-1;
    double loglambda = s.fe + s.R[t] + s.betaX[data.mbrg[u]-1] * s.X(t, data.mbrg[u]-1);
    lr += Util::loglik_pois(data.cases(t,u), data.n[u], loglambda+s.U(u,p), loglambda+prop);
  }
  return lr;
//  tps <- length(R)
//  lambda_curr <- n[j]*exp(fe+R+curr+rep(betaX[mbrg[j]],tps)*X[,mbrg[j]])
//  lambda_prop <- lambda_curr * exp(prop-curr)
//  sum(cases[,j] * (prop - curr) - lambda_prop + lambda_curr)
}

double sum_u_squared(const NumericMatrix &U, const NumericMatrix &nb) {
  double sum = 0;
  for (int p = 0; p < U.ncol(); p++) {
    for (int i = 0; i < U.nrow(); i++) {
      for (int j = 1; j < nb(i,0)+1; j++) {
        double d = U(i,p) - U(nb(i, j)-1,p);
        sum += d*d;
      }
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
  s.kU = Util::rgamma(aU + 0.5*(s.U.nrow()-1)*s.U.ncol(), bU + 0.5*sum_u_squared(s.U, data.nb));

  for (int p = 0; p < s.U.ncol(); p++) {
    const NumericVector &p2t = data.p2t[p];
    for (int u = 0; u < s.U.nrow(); u++) {
      double proposal, ap;
      if (i % 2 == 0) { // propose from the prior
        double sumU = 0;
        for (int j = 1; j < data.nb(u,0) + 1; j++)
          sumU += s.U(data.nb(u,j)-1, p);

        proposal = R::rnorm(sumU/data.nb(u,0), 1/::sqrt(s.kU * data.nb(u,0)));
        ap = u_likelihood(data, s, proposal, u, p);
      } else { // M-H random walk proposal
        proposal = R::rnorm(s.U(u,p), sigmaU);
        double sumU = 0;
        for (int j = 1; j < data.nb(u,0) + 1; j++) {
          double d1 = s.U(data.nb(u,j)-1,p) - proposal;
          double d2 = s.U(data.nb(u,j)-1,p) - s.U(u,p);
          sumU += d1*d1 - d2*d2;
        }
        ap = u_likelihood(data, s, proposal, u, p) - s.kU * sumU / 2;
      }

      // multiply this by the proposal ratio, as we need to adjust R accordingly.
      // If you think about the R precision matrix then that involves the two
      // timepoints either side of the period division(s).
      double r_prior_ratio = 0;
      double delta_U = (proposal - s.U(u,p)) / s.U.nrow();
      if (p > 0) {
        int t = p2t[0]-1;
        r_prior_ratio += s.kR * delta_U * (s.R[t-2] - 3*s.R[t-1] + 3*s.R[t] - s.R[t+1] + delta_U);
      }
      if (p < s.U.ncol()-1) {
        int t = p2t[p2t.length()-1]-1;
        r_prior_ratio += s.kR * delta_U * (s.R[t+2] - 3*s.R[t+1] + 3*s.R[t] - s.R[t-1] + delta_U);
      }
      ap -= r_prior_ratio;

      double un = R::unif_rand();
      if (ap >= 0 || un <= ::exp(ap)) {
        s.U(u,p) = proposal;
        // also update R and U to keep constant mean
        for (int j = 0; j < p2t.length(); j++) {
          int t = p2t[j]-1;
          s.R[t] += delta_U;
        }
        s.U(_,p) = s.U(_,p) - delta_U;
        s.acceptU[i%2]++;
      } else {
        s.rejectU[i%2]++;
      }
    }
  }
}
