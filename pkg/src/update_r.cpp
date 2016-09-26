#include "update.h"
#include "util.h"

using namespace Rcpp;

double sum_r_squared(const NumericVector &R) {
  double sum = 0;
  for (int i = 0; i < R.length()-2; i++) {
    double r = R[i] - 2*R[i+1] + R[i+2];
    sum += r * r;
  }
  return sum;
}

double r_likelihood(const Data &data,
                    const State &s,
                         NumericVector prop,
                         int j) {
  const int n_u = s.U.length();
  double lr = 0;
  for (int u = 0; u < n_u; u++) {
    double l_u = s.fe + s.U[u];
    for (int i = 0; i < prop.length(); i++) {
      double loglambda = l_u + s.X(i+j, data.mbrg[u]-1)*s.betaX[data.mbrg[u]-1];
      lr += data.cases(i+j,u) * (prop[i] - s.R[i+j])
        - data.n[u] * (::exp(loglambda + prop[i]) - ::exp(loglambda + s.R[i+j]));
    }
  }
  return lr;
}

inline double square_diff(double a, double b) {
  return a*a - b*b;
}

void update_r(const Data &data,
                int i,
                State &s,
                List prior,
                List control) {

  double aR = prior["aR"];
  double bR = prior["bR"];
  double sigmaR = control["sigmaR"];
  List Rbefore  = control["Rbefore"];
  List Rafter   = control["Rafter"];
  List Rsigma   = control["Rsigma"];
  int max_methods = Rbefore.size() / 5;

  // Gibbs update for kR
  s.kR = Util::rgamma(aR + 0.5*(s.R.length()-2), bR+0.5*sum_r_squared(s.R));

  int method = i % (1+max_methods);
  int endmethod = Util::rbernoulli(0.5);

  int j = 0; // start of update block
  while (j < s.R.length()) {
    NumericVector proposal;
    double ap;
    if (method >= max_methods || (endmethod == 0 && (j < 2 || j > s.R.length()-3))) {
      // Metropolis Hastings proposal step to update R.
      double prop = R::rnorm(s.R[j], sigmaR);
      double prior_ratio = 0;
      if (j > 1)
        prior_ratio += 0.5 * s.kR * square_diff(s.R[j-2]-2*s.R[j-1]+s.R[j], s.R[j-2]-2*s.R[j-1]+prop);
      if (j > 0 && j < s.R.length()-1)
        prior_ratio += 0.5 * s.kR * square_diff(s.R[j-1]-2*s.R[j]+s.R[j+1], s.R[j-1]-2*prop+s.R[j+1]);
      if (j < s.R.length()-2)
        prior_ratio += 0.5 * s.kR * square_diff(s.R[j]-2*s.R[j+1]+s.R[j+2], prop-2*s.R[j+1]+s.R[j+2]);
      proposal = prop;
      ap = r_likelihood(data, s, proposal, j) + prior_ratio;
    } else {
      // Conditional Prior Proposal step to update R
/*      if (j == 0)
        proposal = R::rnorm(2*s.R[1]-s.R[2], ::sqrt(1/s.kR));
      else if (j == 1)
        proposal = R::rnorm(0.4*s.R[0]+0.8*s.R[2]-0.2*s.R[3], ::sqrt(0.2/s.kR));
      else if (j == s.R.length() - 2) {
        // this case, and the above are derived from R[j] = 2R[j-1] - R[j-2]
        // combined with twice (for some reason)    2R[j] = R[j+1] + R[j-1]
        // resulting in                             5R[j] = 4R[j-1] + 2R[j+1] - R[j-2]
        proposal = R::rnorm(-0.2*s.R[j-2]+0.8*s.R[j-1]+0.4*s.R[j+1], ::sqrt(0.2/s.kR));
      }
      else if (j == s.R.length() - 1)
        proposal = R::rnorm(-s.R[j-2]+2*s.R[j-1], ::sqrt(1/s.kR));
      else */{
        // first check if we're going to hit the end
        int size = as<NumericMatrix>(Rbefore[method*5]).nrow();
        if (j + size > s.R.length()) { // hit the end, so shuffle down a bit
          j = s.R.length() - size;
        }
        // now work out the appropriate method. We use 2 unless we're at the ends
        int o = 2; // default method for j > 2
        if (j < 2) o = j;
        if (j + size > s.R.length() - 2) o = 4 - (s.R.length() - (j + size));
        NumericMatrix rbe = Rbefore[method*5 + o];
        NumericMatrix raf = Rafter[method*5 + o];
        NumericMatrix rsigma = Rsigma[method*5 + o];
        NumericVector mu(rbe.nrow());
        for (int i = 0; i < mu.length(); i++) {
          for (int l = 0; l < rbe.ncol(); l++) {
            mu[i] += rbe(i,l) * s.R[j - rbe.ncol() + l];
          }
          for (int l = 0; l < raf.ncol(); l++) {
            mu[i] += raf(i,l) * s.R[j + raf.nrow() + l];
          }
        }
        proposal = Util::rmvnorm(mu, rsigma/::sqrt(s.kR));
      }
      ap = r_likelihood(data, s, proposal, j);
    }
    double un = R::unif_rand();
    if (ap >= 0 | un <= ::exp(ap)) {
      std::copy(proposal.begin(), proposal.end(), s.R.begin() + j);
      s.acceptR[method]++;
    } else {
      s.rejectR[method]++;
    }
    j += proposal.length();
  }

  double meanR = mean(s.R);
  s.fe += meanR;
  s.R = s.R - meanR;
}
