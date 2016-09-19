#include <Rcpp.h>
#include <cmath>
#include "util.h"

using namespace Rcpp;

// [[Rcpp::export]]
double sum_r_squared(NumericVector R) {
  double sum = 0;
  for (int i = 0; i < R.length()-2; i++) {
    double r = R[i] - 2*R[i+1] + R[i+2];
    sum += r * r;
  }
  return sum;
}

double r_likelihood(NumericMatrix cases,
                         NumericVector n,
                         NumericVector mbrg,
                         double fe,
                         NumericVector R,
                         NumericVector U,
                         NumericMatrix X,
                         NumericVector betaX,
                         NumericVector prop,
                         int j) {
  const int n_u = U.length();
  double lr = 0;
  for (int u = 0; u < n_u; u++) {
    double l_u = fe + U[u];
    for (int i = 0; i < prop.length(); i++) {
      double loglambda = l_u + X(i+j, mbrg[u]-1)*betaX[mbrg[u]-1];
      lr += cases(i+j,u) * (prop[i] - R[i+j])
        - n[u] * (::exp(loglambda + prop[i]) - ::exp(loglambda + R[i+j]));
    }
  }
  return lr;
}

inline double square_diff(double a, double b) {
  return a*a - b*b;
}

// [[Rcpp::export]]
Rcpp::List update_r(NumericMatrix cases,
                NumericVector n,
                NumericVector mbrg,
                int i,
                List state,
                List prior,
                List Rmu,
                List Rsigma) {

  RNGScope scope;
  List out = state;

  // extract information from our state
  NumericVector R = state["R"];
  double       fe = state["fe"];
  NumericVector U = state["U"];
  NumericMatrix X = state["X"];
  NumericVector betaX = state["betaX"];
  NumericVector acceptR = state["acceptR"];
  NumericVector rejectR = state["rejectR"];

  double aR = prior["aR"];
  double bR = prior["bR"];
  double sigmaR = prior["sigmaR"];

  // Gibbs update for kR
  double kR = Util::rgamma(aR + 0.5*(R.length()-2), bR+0.5*sum_r_squared(R));

  int method = i % (1+Rmu.size());
  int endmethod = Util::rbernoulli(0.5);

  int j = 0; // start of update block
  while (j < R.length()) {
    NumericVector proposal;
    double ap;
    if (method >= Rmu.size() || (endmethod == 0 && (j < 2 || j > R.length()-3))) {
      // Metropolis Hastings proposal step to update R.
      double prop = R::rnorm(R[j], sigmaR);
      double prior_ratio = 0;
      if (j > 1)
        prior_ratio += 0.5 * kR * square_diff(R[j-2]-2*R[j-1]+R[j], R[j-2]-2*R[j-1]+prop);
      if (j > 0 && j < R.length()-1)
        prior_ratio += 0.5 * kR * square_diff(R[j-1]-2*R[j]+R[j+1], R[j-1]-2*prop+R[j+1]);
      if (j < R.length()-2)
        prior_ratio += 0.5 * kR * square_diff(R[j]-2*R[j+1]+R[j+2], prop-2*R[j+1]+R[j+2]);
      proposal = prop;
      ap = r_likelihood(cases, n, mbrg, fe, R, U, X, betaX, proposal, j) + prior_ratio;
    } else {
      // Conditional Prior Proposal step to update R
      if (j == 0)
        proposal = R::rnorm(2*R[1]-R[2], ::sqrt(1/kR));
      else if (j == 1)
        proposal = R::rnorm(0.4*R[0]+0.8*R[2]-0.2*R[3], ::sqrt(0.2/kR));
      else if (j == R.length() - 2) // TODO: I think the +0.2 should be +0.4...
        proposal = R::rnorm(-0.2*R[j-2]+0.8*R[j-1]+0.2*R[j+1], ::sqrt(0.2/kR));
      else if (j == R.length() - 1)
        proposal = R::rnorm(-R[j-2]+2*R[j-1], ::sqrt(1/kR));
      else {
        NumericMatrix rmu = Rmu[method];
        NumericMatrix rsigma = Rsigma[method];
        if (j + rmu.nrow() > R.length() - 2) {
          j = R.length() - 2 - rmu.nrow();
        }
        NumericVector mu = no_init(rmu.nrow());
        for (int i = 0; i < mu.length(); i++)
          mu[i] = rmu(i,0)*R[j-2]+rmu(i,1)*R[j-1]+rmu(i,2)*R[j+rmu.nrow()]+rmu(i,3)*R[j+rmu.nrow()+1];
        proposal = Util::rmvnorm(mu, rsigma/::sqrt(kR));
      }
      ap = r_likelihood(cases, n, mbrg, fe, R, U, X, betaX, proposal, j);
    }
    double un = R::unif_rand();
    if (ap >= 0 | un <= ::exp(ap)) {
      std::copy(proposal.begin(), proposal.end(), R.begin() + j);
      acceptR[method]++;
    } else {
      rejectR[method]++;
    }
    j += proposal.length();
  }
  out["R"] = R;
  out["kR"] = kR;

  out["acceptR"] = acceptR;
  out["rejectR"] = rejectR;
  return out;
}
