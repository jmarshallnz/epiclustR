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
        loglr += Util::loglik_pois(data.cases(t,u), data.n[u], loglambda0, loglambda1);
      }
      double pr = s.pX / (::exp(-loglr) * (1-s.pX) + s.pX);
      s.X(t, r) = Util::rbernoulli(pr);
    }
  }
//  squashProd(dpois(cases,rep(n,each=nrow(X))*exp(fe+rep(R,ncol(n))+rep(U,each=nrow(X))+0*rep(betaX[mbrg],each=nrow(X))), log=TRUE)-
//             dpois(cases,rep(n,each=nrow(X))*exp(fe+rep(R,ncol(n))+rep(U,each=nrow(X))+1*rep(betaX[mbrg],each=nrow(X))), log=TRUE))

//  lenR <- length(R)
//  lambda0 <- rep(n,each=lenR)*exp(fe+rep(R,ncol(n))+rep(U,each=lenR))
//  lambda1 <- lambda0 * exp(rep(betaX[mbrg],each=lenR))
//  squashProd(cases * (0 - rep(betaX[mbrg],each=lenR)) - lambda0 + lambda1)

//  loglambda0 <- fe+rep(R,ncol(n))+rep(U,each=lenR)
//  loglambda1 <- lambda0 + rep(betaX[mbrg],each=lenR)
//  squashProd(n * (exp(loglambda1) - exp(loglambda1)) - cases * (loglambda1 - loglambda0))
}

void sample_x2(const Data &data, State &s) {
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
        if (t > 0)
          loglambda0 += s.X(t-1,r) * s.betaX[r];
        double loglambda1 = loglambda0 + s.betaX[r];
        loglr += Util::loglik_pois(data.cases(t,u), data.n[u], loglambda0, loglambda1);
        if (t+1 < n_t) {
          int p2 = data.t2p[t+1]-1;
          loglambda0 = s.fe + s.R[t+1] + s.U(u,p2) + s.X(t+1,r) * s.betaX[r];
          loglambda1 = loglambda0 + s.betaX[r];
          loglr += Util::loglik_pois(data.cases(t+1,u), data.n[u], loglambda0, loglambda1);
        }
      }
      double pr = s.pX / (::exp(-loglr) * (1-s.pX) + s.pX);
      s.X(t, r) = Util::rbernoulli(pr);
    }
  }
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
      lr += Util::loglik_pois(data.cases(t,u), data.n[u], loglambda0, loglambda1);
    }
  }
  return lr;
  // Xj <- X[rep(j-1,each=tps)*tps+rep(1:tps,lwch[j])]
  // lambda_curr <- rep(n[wch[[j]]],each=tps)*exp(fe+rep(R,lwch[j])+rep(U[wch[[j]]],each=tps)+Xj*curr)
  // lambda_prop <- lambda_curr * exp(Xj*(prop-curr))
  // sum(cases[,wch[[j]]] * X * (prop - curr) - lambda_prop + lambda_curr)
}

double betax_likelihood2(const Data &data,
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
      if (t > 0) {
        loglambda0 += s.X(t-1,j)*s.betaX[j];
        loglambda1 += s.X(t-1,j)*prop;
      }
      lr += Util::loglik_pois(data.cases(t,u), data.n[u], loglambda0, loglambda1);
    }
  }
  return lr;
  // Xj <- X[rep(j-1,each=tps)*tps+rep(1:tps,lwch[j])]
  // lambda_curr <- rep(n[wch[[j]]],each=tps)*exp(fe+rep(R,lwch[j])+rep(U[wch[[j]]],each=tps)+Xj*curr)
  // lambda_prop <- lambda_curr * exp(Xj*(prop-curr))
  // sum(cases[,wch[[j]]] * X * (prop - curr) - lambda_prop + lambda_curr)

  // lambda_curr0 <- n[wch[[j]]] * exp(fe+rep(R[1],lwch[j])+U[wch[[j]]]+X[1,j]*curr)
  // lambda_prop0 <- n[wch[[j]]] * exp(fe+rep(R[1],lwch[j])+U[wch[[j]]]+X[1,j]*prop)
  // lr += dpois(cases(1, wch[[j]]), lambda_prop0)-
  //       dpois(cases(1, wch[[j]]), lambda_curr0)
  // lambda_curr1 <- rep(n[wch[[j]]],each=tps-1) * exp(fe+rep(R[2:tps],lwch[j])+rep(U[wch[[j]]],each=tps-1)+(X[rep((j-1)*tps,each=tps-1)+rep(2:tps,lwch[j])]+X[rep((j-1)*tps,each=tps-1)+rep(1:(tps-1),lwch[j])])*curr)
  // lambda_prop1 <- rep(n[wch[[j]]],each=tps-1) * exp(fe+rep(R[2:tps],lwch[j])+rep(U[wch[[j]]],each=tps-1)+(X[rep((j-1)*tps,each=tps-1)+rep(2:tps,lwch[j])]+X[rep((j-1)*tps,each=tps-1)+rep(1:(tps-1),lwch[j])])*prop)
  // lr += dpois(cases(2:tps, wch[[j]]), lambda_prop1)-
  //       dpois(cases(2:tps, wch[[j]]), lambda_curr1)
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
