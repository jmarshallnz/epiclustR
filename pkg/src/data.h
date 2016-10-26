#pragma once

#include <Rcpp.h>

class Data {
  public:
    Rcpp::NumericMatrix cases;
    Rcpp::NumericVector n;
    Rcpp::NumericVector mbrg;
    Rcpp::NumericMatrix nb;
    Rcpp::List          rgmb;

    Data(Rcpp::List d) :
      cases(Rcpp::as<Rcpp::NumericMatrix>(d["cases"])),
      n(Rcpp::as<Rcpp::NumericVector>(d["popn"])),
      mbrg(Rcpp::as<Rcpp::NumericVector>(d["mbrg"])),
      nb(Rcpp::as<Rcpp::NumericMatrix>(d["nb"])),
      rgmb(Rcpp::as<Rcpp::List>(d["rgmb"])) {
    }
};
