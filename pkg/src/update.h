#pragma once

#include <Rcpp.h>
#include "data.h"

Rcpp::List update_r(const Data &data,
                    int i,
                    Rcpp::List state,
                    Rcpp::List prior,
                    Rcpp::List control);

Rcpp::List update_u(const Data &data,
                    int i,
                    Rcpp::List state,
                    Rcpp::List prior,
                    Rcpp::List control);

Rcpp::List update_x(const Data &data,
                    Rcpp::List state,
                    Rcpp::List prior,
                    Rcpp::List control);
