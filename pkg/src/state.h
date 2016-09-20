#pragma once

#include <Rcpp.h>

class State {
  public:
    double              fe;
    Rcpp::NumericVector R;
    Rcpp::NumericVector U;
    Rcpp::IntegerMatrix X;
    Rcpp::NumericVector betaX;
    double              kR;
    double              kU;
    double              pX;
    Rcpp::NumericVector acceptR;
    Rcpp::NumericVector rejectR;
    Rcpp::NumericVector acceptU;
    Rcpp::NumericVector rejectU;
    double acceptX;
    double rejectX;

    State(Rcpp::List s) :
      fe(Rcpp::as<double>(s["fe"])),
      R(Rcpp::as<Rcpp::NumericVector>(s["R"])),
      U(Rcpp::as<Rcpp::NumericVector>(s["U"])),
      X(Rcpp::as<Rcpp::IntegerMatrix>(s["X"])),
      betaX(Rcpp::as<Rcpp::NumericVector>(s["betaX"])),
      kR(Rcpp::as<double>(s["kR"])),
      kU(Rcpp::as<double>(s["kU"])),
      pX(Rcpp::as<double>(s["pX"])),
      acceptR(Rcpp::as<Rcpp::NumericVector>(s["acceptR"])),
      rejectR(Rcpp::as<Rcpp::NumericVector>(s["rejectR"])),
      acceptU(Rcpp::as<Rcpp::NumericVector>(s["acceptU"])),
      rejectU(Rcpp::as<Rcpp::NumericVector>(s["rejectU"])),
      acceptX(Rcpp::as<double>(s["acceptX"])),
      rejectX(Rcpp::as<double>(s["rejectX"])) {
    }

    Rcpp::List toList(Rcpp::List i) const {
      Rcpp::List s = i;
      s["fe"] = fe;
      s["R"] = R;
      s["U"] = U;
      s["X"] = X;
      s["betaX"] = betaX;
      s["kR"] = kR;
      s["kU"] = kU;
      s["pX"] = pX;
      s["acceptR"] = acceptR;
      s["rejectR"] = rejectR;
      s["acceptU"] = acceptU;
      s["rejectU"] = rejectU;
      s["acceptX"] = acceptX;
      s["rejectX"] = rejectX;
      return s;
    }
};
