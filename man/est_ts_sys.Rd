% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/two_phase_estimators.R
\name{est_ts_sys}
\alias{est_ts_sys}
\title{The population is separated into a set of contiguous PSUs (clusters/strips),
which are sampled in the first stage. SSUs are selected systematically across
the selected PSUs in the second stage. Adapted to systematic sampling.}
\usage{
est_ts_sys(dt_samp, id_psu = "id", y = quote(y), a_psu, a_ssu, NI, nI,
  N_i = quote(N_i))
}
\arguments{
\item{id_psu}{Name of the variable identifying the PSUs. Must be present in
both data.tables.}

\item{y}{Quoted name of the variable storing the observed values of the study
variable.}

\item{a_psu}{PSU intervall, every a-th PSU is selected into the sample.}

\item{a_ssu}{SSU intervall, every a-th SSU is selected into the sample.}

\item{NI}{Number of PSUs}

\item{nI}{PSU sample size}

\item{N_i}{Optional. PSU size in case PSUs are of different size. Needed for
the successive differences variance estimator.}

\item{dt_sample}{A data.table object containing information about the PSUs and
SSUs selected in the two sampling stages.}
}
\value{
A named list containing estimates for the total and its variance
 estimate
}
\description{
The population is separated into a set of contiguous PSUs (clusters/strips),
which are sampled in the first stage. SSUs are selected systematically across
the selected PSUs in the second stage. Adapted to systematic sampling.
}
\details{
See Särndal et al. (1992), equations 4.3.21 for the estimation of
 the population total, slightly adapted to case of systematic sampling. See
 Ene et al. (2016, equation 13) for variance estimation using successive
 differences.
}
\references{
Ene, L. T., E. Næsset and T. Gobakken (2016). "Simulation-based assessment
 of sampling strategies for large-area biomass estimation using wall-to-wall
 and partial coverage airborne laser scanning surveys." Remote Sensing of
 Environment 176: 328-340.

 Särndal, C. E., B. Swensson and J. Wretman (1992). Model assisted survey
 sampling. New York, Springer-Verlag.
}

