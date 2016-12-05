#pragma once

#include <Rcpp.h>

class Data {
  public:
    Rcpp::NumericMatrix cases;
    Rcpp::NumericMatrix n;
    Rcpp::NumericVector mbrg;
    Rcpp::NumericMatrix nb;
    Rcpp::List          rgmb;
    Rcpp::NumericVector t2p;
    Rcpp::List          p2t;

    Data(Rcpp::List d) :
      cases(Rcpp::as<Rcpp::NumericMatrix>(d["cases"])),
      n(Rcpp::as<Rcpp::NumericMatrix>(d["popn"])),
      mbrg(Rcpp::as<Rcpp::NumericVector>(d["mbrg"])),
      nb(Rcpp::as<Rcpp::NumericMatrix>(d["nb"])),
      rgmb(Rcpp::as<Rcpp::List>(d["rgmb"])),
      t2p(Rcpp::as<Rcpp::NumericVector>(d["t2p"])),
      p2t(Rcpp::as<Rcpp::List>(d["p2t"])) {
    }
};
