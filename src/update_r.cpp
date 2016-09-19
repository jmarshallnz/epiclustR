#include <Rcpp.h>
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
        - n[u] * (exp(loglambda + prop[i]) - exp(loglambda + R[i+j]));
    }
  }
  return lr;
}

double r_likelihood(NumericMatrix cases,
                    NumericVector n,
                    NumericVector mbrg,
                    double fe,
                    NumericVector U,
                    NumericMatrix X,
                    NumericVector betaX,
                    double curr,
                    double prop,
                    int j) {
  const int n_u = U.length();
  double lr = 0;
  for (int u = 0; u < n_u; u++) {
    double loglambda = fe + U[u] + X(j, mbrg[u]-1)*betaX[mbrg[u]-1];
    lr += cases(j,u) * (prop - curr)
      - n[u] * (exp(loglambda + prop) - exp(loglambda + curr));
  }
  return lr;
}

inline double square_diff(double a, double b) {
  return a*a - b*b;
}

// [[Rcpp::export]]
Rcpp::List update_r_mh(NumericMatrix cases,
                    NumericVector n,
                    NumericVector mbrg,
                    List state,
                    List prior,
                    int j) {

  RNGScope scope;
  Rcpp::List out;

  // extract information from our state
  NumericVector R = state["R"];
  double       kR = state["kR"];
  double       fe = state["fe"];
  NumericVector U = state["U"];
  NumericMatrix X = state["X"];
  NumericVector betaX = state["betaX"];

  double sigmaR = prior["sigmaR"];

  double proposal = R::rnorm(R[j], sigmaR);

  double ap = r_likelihood(cases, n, mbrg, fe, U, X, betaX, R[j], proposal, j);
  double prior_ratio = 0;
  if (j > 1)
    prior_ratio += 0.5 * kR * square_diff(R[j-2]-2*R[j-1]+R[j], R[j-2]-2*R[j-1]+proposal);
  if (j > 0 && j < R.length()-1)
    prior_ratio += 0.5 * kR * square_diff(R[j-1]-2*R[j]+R[j+1], R[j-1]-2*proposal+R[j+1]);
  if (j < R.length()-2)
    prior_ratio += 0.5 * kR * square_diff(R[j]-2*R[j+1]+R[j+2], proposal-2*R[j+1]+R[j+2]);

  out["proposal"] = proposal;
  out["ap"] = ap + prior_ratio;

  return out;
}

// [[Rcpp::export]]
Rcpp::List update_r_cond(NumericMatrix cases,
                       NumericVector n,
                       NumericVector mbrg,
                       List state,
                       NumericMatrix Rmu,
                       NumericMatrix Rsigma,
                       int j) {

  RNGScope scope;
  Rcpp::List out;

  // extract information from our state
  NumericVector R = state["R"];
  double       kR = state["kR"];
  double       fe = state["fe"];
  NumericVector U = state["U"];
  NumericMatrix X = state["X"];
  NumericVector betaX = state["betaX"];

  NumericVector proposal;
  if (j == 0)
    proposal = R::rnorm(2*R[1]-R[2], sqrt(1/kR));
  else if (j == 1)
    proposal = R::rnorm(0.4*R[0]+0.8*R[2]-0.2*R[3], sqrt(0.2/kR));
  else if (j == R.length() - 2) // TODO: I think the +0.2 should be +0.4...
    proposal = R::rnorm(-0.2*R[j-2]+0.8*R[j-1]+0.2*R[j+1], sqrt(0.2/kR));
  else if (j == R.length() - 1)
    proposal = R::rnorm(-R[j-2]+2*R[j-1], sqrt(1/kR));
  else {
    int k = j + Rmu.nrow() - 1;
    if (k > R.length()-3) {
      k = R.length() - 3;
      j = k - Rmu.nrow() + 1;
    }
    NumericVector mu = no_init(Rmu.nrow());
    for (int i = 0; i < mu.length(); i++)
      mu[i] = Rmu(i,0)*R[j-2]+Rmu(i,1)*R[j-1]+Rmu(i,2)*R[k+1]+Rmu(i,3)*R[k+2];
    proposal = Util::rmvnorm(mu, Rsigma/sqrt(kR));
  }
  double ap = r_likelihood(cases, n, mbrg, fe, R, U, X, betaX, proposal, j);

  out["proposal"] = proposal;
  out["ap"] = ap;

  return out;
}

// [[Rcpp::export]]
Rcpp::List update_r(NumericMatrix cases,
                NumericVector n,
                int i,
                List state,
                List prior) {

  RNGScope scope;

  // extract information from our state
  NumericVector R = state["R"];

  double aR = prior["aR"];
  double bR = prior["bR"];

  // Gibbs update for kR
  double kR = Util::rgamma(aR + 0.5*(R.length()-2), bR+0.5*sum_r_squared(R));
  state["kR"] = kR;

  return state;
}
