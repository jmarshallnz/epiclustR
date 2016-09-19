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

Rcpp::IntegerMatrix sample_x(NumericMatrix cases,
                                  NumericVector n,
                                  double fe,
                                  NumericVector R,
                                  NumericVector U,
                                  NumericVector betaX,
                                  double pX,
                                  List rgmb) {
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

double sample_px(IntegerMatrix X, double aX, double bX) {
  double sumX = sum(X);
  return R::rbeta(aX+sumX, bX+X.nrow()*X.ncol()-sumX);
}

double betax_likelihood(NumericMatrix cases,
                             NumericVector n,
                             double fe,
                             NumericVector R,
                             NumericVector U,
                             IntegerMatrix X,
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
      lr += cases(t,u) * X(t,j) * (prop - curr)
        - n[u] * (exp(loglambda + X(t,j)*prop) - exp(loglambda + X(t,j)*curr));
    }
  }
  return lr;
  // Xj <- X[rep(j-1,each=tps)*tps+rep(1:tps,lwch[j])]
  // lambda_curr <- rep(n[wch[[j]]],each=tps)*exp(fe+rep(R,lwch[j])+rep(U[wch[[j]]],each=tps)+Xj*curr)
  // lambda_prop <- lambda_curr * exp(Xj*(prop-curr))
  // sum(cases[,wch[[j]]] * X * (prop - curr) - lambda_prop + lambda_curr)
}

// [[Rcpp::export]]
Rcpp::List update_x(NumericMatrix cases,
                    NumericVector n,
                    List rgmb,
                    List state,
                    List prior) {

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
  IntegerMatrix X = sample_x(cases, n, fe, R, U, betaX, pX, rgmb);
  out["X"]  = X;

  // sample pX
  pX = sample_px(X, prior["aX"], prior["bX"]);
  out["pX"] = pX;

  // sample betaX
  double sigmaX = prior["sigmaX"];
  double abetaX = prior["abetaX"];
  double bbetaX = prior["bbetaX"];

  for (int r = 0; r < betaX.length(); r++) {
    double proposal = R::rnorm(betaX[r], sigmaX);
    if (proposal <= 0) {
       rejectX++;
    } else {
      double prior_ratio = (abetaX - 1) * (log(proposal) - log(betaX[r])) - (proposal - betaX[r])*bbetaX;
      double ap = betax_likelihood(cases, n, fe, R, U, X, rgmb[r], betaX[r], proposal, r) + prior_ratio;
      double un = R::unif_rand();
      if (ap >= 0 || un <= exp(ap)) {
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
