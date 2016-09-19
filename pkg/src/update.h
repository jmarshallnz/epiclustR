#pragma once

#include <Rcpp.h>

Rcpp::List update_r(Rcpp::NumericMatrix cases,
                    Rcpp::NumericVector n,
                    Rcpp::NumericVector mbrg,
                    int i,
                    Rcpp::List state,
                    Rcpp::List prior,
                    Rcpp::List control);

Rcpp::List update_u(Rcpp::NumericMatrix cases,
                    Rcpp::NumericVector n,
                    Rcpp::NumericVector mbrg,
                    Rcpp::NumericMatrix nb,
                    int i,
                    Rcpp::List state,
                    Rcpp::List prior,
                    Rcpp::List control);

Rcpp::List update_x(Rcpp::NumericMatrix cases,
                    Rcpp::NumericVector n,
                    Rcpp::List rgmb,
                    Rcpp::List state,
                    Rcpp::List prior,
                    Rcpp::List control);
