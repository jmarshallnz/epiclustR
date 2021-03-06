% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostics.R
\name{rhat}
\alias{rhat}
\title{Compute the split-R-hat (a measure of convergence) for a variable from posterior samples.}
\usage{
rhat(x)
}
\arguments{
\item{x}{A matrix or vector of posterior output for a variable. Rows are samples, columns are chains.}
}
\value{
The split-R-hat for the variable in question. Should be close to 1 on convergence.
}
\description{
Computation taken from the `rstan` package. Included here so that the package
can be as self-contained as possible.
}
