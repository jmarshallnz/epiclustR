% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{fit_model}
\alias{fit_model}
\title{Fits the epiclustR model for the given dataset}
\usage{
fit_model(data, prior, control, seed = NULL)
}
\arguments{
\item{data}{the data for the model}

\item{prior}{details on the priors for the model}

\item{control}{details on the control for the model}

\item{seed}{the random number seed to set. Defaults to NULL (a random seed)}
}
\value{
the posterior for this chain
}
\description{
Fits the epiclustR model for the given dataset
}
