#include <Rcpp.h>
#include "update.h"
#include "data.h"

using namespace Rcpp;

//' Update the MCMC simulation
//'
//' @param data a list containing the data
//' @param state the current state of the Markov chain
//' @param prior the priors used for the model
//' @param control the control parameters for the MCMC algorithm (e.g. tuning, blocks)
//' @return a list containing the updated state
//' @export
// [[Rcpp::export]]
Rcpp::List update(List data,
                  List state,
                  List prior,
                  List control) {

  // setup our objects
  Data d(data);
  State s(state);

  int thinning = control["thinning"];

  for (int i = 0; i < thinning; i++) {
    update_r(d, i, s, prior, control);
    update_u(d, i, s, prior, control);
    update_x(d, s, prior, control);
  }

  return s.toList(state);
}
