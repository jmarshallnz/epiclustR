#pragma once

#include <Rcpp.h>
#include "data.h"
#include "state.h"

void update_r(const Data &data,
                    int i,
                    State &state,
                    Rcpp::List prior,
                    Rcpp::List control);

void update_u(const Data &data,
                    int i,
                    State &state,
                    Rcpp::List prior,
                    Rcpp::List control);

void update_x(const Data &data,
                    State &state,
                    Rcpp::List prior,
                    Rcpp::List control);
