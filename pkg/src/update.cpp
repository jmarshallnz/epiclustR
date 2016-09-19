#include <Rcpp.h>
#include "update.h"

using namespace Rcpp;

//' Update the MCMC simulation
//'
//' @param data a list containing the data
//' @param i the current iteration
//' @param state the current state of the Markov chain
//' @param prior the priors used for the model
//' @param control the control parameters for the MCMC algorithm (e.g. tuning, blocks)
//' @return a list containing the updated state
//' @export
// [[Rcpp::export]]
Rcpp::List update(List data,
                  int i,
                  List state,
                  List prior,
                  List control) {

  // setup our objects
  NumericMatrix cases = data["cases"];
  NumericVector     n = data["popn"];
  NumericVector  mbrg = data["mbrg"];
  NumericMatrix    nb = data["nb"];
  List           rgmb = data["rgmb"];

  List s = state;
  s = update_r(cases, n, mbrg, i, s, prior, control);
  s = update_u(cases, n, mbrg, nb, i, s, prior, control);
  s = update_x(cases, n, rgmb, s, prior, control);
  return s;
}
