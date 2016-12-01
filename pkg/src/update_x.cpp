#include "update.h"
#include "util.h"

using namespace Rcpp;

void sample_x(const Data &data, State &s) {
  const int n_t = s.R.length();
  const int n_r = data.rgmb.length();
  for (int r = 0; r < n_r; r++) {
    const NumericVector &rgmb_r = data.rgmb[r];
    for (int t = 0; t < n_t; t++) {
      double loglr = 0;
      int p = data.t2p[t]-1;
      for (int j = 0; j < rgmb_r.length(); j++) {
        int u = rgmb_r[j]-1;
        double loglambda0 = s.fe + s.R[t] + s.U(u,p);
        double loglambda1 = loglambda0 + s.betaX[r];
        loglr += Util::loglik_pois(data.cases(t,u), data.n(t,u), loglambda0, loglambda1);
      }
      double pr = s.pX / (::exp(-loglr) * (1-s.pX) + s.pX);
      s.X(t, r) = Util::rbernoulli(pr);
    }
  }
//  lenR <- length(R)
//  lambda0 <- rep(n,each=lenR)*exp(fe+rep(R,ncol(n))+rep(U,each=lenR))
//  lambda1 <- lambda0 * exp(rep(betaX[mbrg],each=lenR))
//  squashProd(cases * (0 - rep(betaX[mbrg],each=lenR)) - lambda0 + lambda1)
}

double sample_px(const IntegerMatrix &X, double aX, double bX) {
  double sumX = sum(X);
  return R::rbeta(aX+sumX, bX+X.nrow()*X.ncol()-sumX);
}

double betax_likelihood(const Data &data,
                        const State &s,
                             double prop,
                             int j) {
  const int n_t = s.R.length();
  NumericVector rgmb_j = data.rgmb[j];
  const int n_r = rgmb_j.length();
  double lr = 0;
  for (int t = 0; t < n_t; t++) {
    for (int r = 0; r < n_r; r++) {
      int u = rgmb_j[r]-1;
      int p = data.t2p[t]-1;
      double loglambda0 = s.fe + s.R[t] + s.U(u,p) + s.X(t,j)*s.betaX[j];
      double loglambda1 = s.fe + s.R[t] + s.U(u,p) + s.X(t,j)*prop;
      lr += Util::loglik_pois(data.cases(t,u), data.n(t,u), loglambda0, loglambda1);
    }
  }
  return lr;
  // Xj <- X[rep(j-1,each=tps)*tps+rep(1:tps,lwch[j])]
  // lambda_curr <- rep(n[wch[[j]]],each=tps)*exp(fe+rep(R,lwch[j])+rep(U[wch[[j]]],each=tps)+Xj*curr)
  // lambda_prop <- lambda_curr * exp(Xj*(prop-curr))
  // sum(cases[,wch[[j]]] * X * (prop - curr) - lambda_prop + lambda_curr)
}

void update_x(const Data &data,
                    State &s,
                    List prior,
                    List control) {

  RNGScope scope;

  // sample X
  sample_x(data, s);

  // sample pX
  s.pX = sample_px(s.X, prior["aX"], prior["bX"]);

  // sample betaX
  double abetaX = prior["abetaX"];
  double bbetaX = prior["bbetaX"];
  double sigmaX = control["sigmaX"];

  for (int r = 0; r < s.betaX.length(); r++) {
    double proposal = R::rnorm(s.betaX[r], sigmaX);
    if (proposal <= 0) {
       s.rejectX++;
    } else {
      double prior_ratio = (abetaX - 1) * (::log(proposal) - ::log(s.betaX[r])) - (proposal - s.betaX[r])*bbetaX;
      double ap = betax_likelihood(data, s, proposal, r) + prior_ratio;
      double un = R::unif_rand();
      if (ap >= 0 || un <= ::exp(ap)) {
        s.betaX[r] = proposal;
        s.acceptX++;
      } else {
        s.rejectX++;
      }
    }
  }
}
