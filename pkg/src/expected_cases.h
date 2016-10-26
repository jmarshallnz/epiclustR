#pragma once

#include <Rcpp.h>
#include "data.h"

Rcpp::NumericMatrix log_case_rate(const Data &data, Rcpp::List state, bool smoothed = false, Rcpp::IntegerVector urange = Rcpp::IntegerVector(0));
