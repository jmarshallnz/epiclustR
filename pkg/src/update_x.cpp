#include "update.h"
#include "util.h"

using namespace Rcpp;

Rcpp::IntegerMatrix sample_x(const Data &data,
                                  double fe,
                                  NumericVector R,
                                  NumericVector U,
                                  NumericVector betaX,
                                  double pX) {
  const int n_t = R.length();
  const int n_r = data.rgmb.length();
  IntegerMatrix X = no_init(n_t, n_r);
  for (int r = 0; r < n_r; r++) {
    for (int t = 0; t < n_t; t++) {
      double loglr = 0;
      NumericVector rgmb_r = data.rgmb[r];
      for (int j = 0; j < rgmb_r.length(); j++) {
        int u = rgmb_r[j]-1;
        double loglambda0 = fe + R[t] + U[u];
        double loglambda1 = loglambda0 + betaX[r];
        loglr += data.n[u] * (::exp(loglambda1) - ::exp(loglambda0)) - data.cases(t,u) * betaX[r];
      }
      double p = pX / (::exp(loglr) * (1-pX) + pX);
      X(t, r) = Util::rbernoulli(p);
    }
  }
  return X;
//  lenR <- length(R)
//  lambda0 <- rep(n,each=lenR)*exp(fe+rep(R,ncol(n))+rep(U,each=lenR))
//  lambda1 <- lambda0 * exp(rep(betaX[mbrg],each=lenR))
//  squashProd(cases * (0 - rep(betaX[mbrg],each=lenR)) - lambda0 + lambda1)
}

double sample_px(IntegerMatrix X, double aX, double bX) {
  double sumX = sum(X);
  return R::rbeta(aX+sumX, bX+X.nrow()*X.ncol()-sumX);
}

double betax_likelihood(const Data &data,
                             double fe,
                             NumericVector R,
                             NumericVector U,
                             IntegerMatrix X,
                             double curr,
                             double prop,
                             int j) {
  const int n_t = R.length();
  NumericVector rgmb_j = data.rgmb[j];
  const int n_r = rgmb_j.length();
  double lr = 0;
  for (int t = 0; t < n_t; t++) {
    for (int r = 0; r < n_r; r++) {
      int u = rgmb_j[r]-1;
      double loglambda = fe + R[t] + U[u];
      lr += data.cases(t,u) * X(t,j) * (prop - curr)
        - data.n[u] * (::exp(loglambda + X(t,j)*prop) - ::exp(loglambda + X(t,j)*curr));
    }
  }
  return lr;
  // Xj <- X[rep(j-1,each=tps)*tps+rep(1:tps,lwch[j])]
  // lambda_curr <- rep(n[wch[[j]]],each=tps)*exp(fe+rep(R,lwch[j])+rep(U[wch[[j]]],each=tps)+Xj*curr)
  // lambda_prop <- lambda_curr * exp(Xj*(prop-curr))
  // sum(cases[,wch[[j]]] * X * (prop - curr) - lambda_prop + lambda_curr)
}

Rcpp::List update_x(const Data &data,
                    List state,
                    List prior,
                    List control) {

  RNGScope scope;

  // grab stuff we need out of the state
  List out = state;
  double fe           = state["fe"];
  NumericVector R     = state["R"];
  NumericVector U     = state["U"];
  NumericVector betaX = state["betaX"];
  double pX           = state["pX"];
  int acceptX         = state["acceptX"];
  int rejectX         = state["rejectX"];

  // sample X
  IntegerMatrix X = sample_x(data, fe, R, U, betaX, pX);
  out["X"]  = X;

  // sample pX
  pX = sample_px(X, prior["aX"], prior["bX"]);
  out["pX"] = pX;

  // sample betaX
  double abetaX = prior["abetaX"];
  double bbetaX = prior["bbetaX"];
  double sigmaX = control["sigmaX"];

  for (int r = 0; r < betaX.length(); r++) {
    double proposal = R::rnorm(betaX[r], sigmaX);
    if (proposal <= 0) {
       rejectX++;
    } else {
      double prior_ratio = (abetaX - 1) * (::log(proposal) - ::log(betaX[r])) - (proposal - betaX[r])*bbetaX;
      double ap = betax_likelihood(data, fe, R, U, X, betaX[r], proposal, r) + prior_ratio;
      double un = R::unif_rand();
      if (ap >= 0 || un <= ::exp(ap)) {
        betaX[r] = proposal;
        acceptX++;
      } else {
        rejectX++;
      }
    }
  }
  out["betaX"]   = betaX;
  out["acceptX"] = acceptX;
  out["rejectX"] = rejectX;

  return out;
}
