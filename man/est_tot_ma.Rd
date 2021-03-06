% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/one_phase_estimators.R
\name{est_tot_ma}
\alias{est_tot_ma}
\title{Estimate the population total using the generalized regression estimator
(GREG)}
\usage{
est_tot_ma(t_pred, res, prob, Nd = NULL)
}
\arguments{
\item{t_pred}{Population total of model predicitions}

\item{res}{Residuals between observed and predicted values in the sample}

\item{prob}{Inclusion probabilities of the design}

\item{Nd}{Domain size. If provided a more precise estimator for the domain
total is used.}
}
\value{
A named vector of three elements containing the estimated total, the
 sample size, and the number of non-zero sample units
}
\description{
This function estimates the population total of a study variable using the
GREG estimator. Domain estimation is possible for the two cases that domain
size is either known or unknown. If known, the appropriate variable needs to
be specified.
}
\references{
See Särndal et al. (1992), eq. 10.5.4 and 10.5.5
}

