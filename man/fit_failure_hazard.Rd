% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nuisance-estimators.R
\name{fit_failure_hazard}
\alias{fit_failure_hazard}
\title{Conditional Failure Hazard Estimator}
\usage{
fit_failure_hazard(
  train_data,
  valid_data,
  learners,
  confounders,
  exposure,
  times
)
}
\arguments{
\item{train_data}{A long \code{data.table} containing the observed data.
\code{train_data} is formatted by \code{\link{unihtee}()} and
\code{\link{tte_data_melt}()}.}

\item{valid_data}{An optional \code{data.table} representing a holdout
dataset of the observed data. It is only used for cross-fitting purposes.
Defaults to \code{NULL}.}

\item{learners}{A \code{\link[sl3]{Stack}}, or other learner class
(inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
learners from \pkg{sl3} to estimate the propensity score model.}

\item{confounders}{A \code{character} vector of column names corresponding to
baseline covariates.}

\item{exposure}{A \code{character} corresponding to the exposure variable.}

\item{times}{A \code{character} giving the name corresponding to the times.}
}
\value{
A named \code{list} of three elements. (1) \code{"estimates"}, the
expected failure hazards for each observation at each time point in
\code{valid_data}, if specified, or \code{train_data} otherwise. (2)
\code{"exp_estimates"}, the expected conditional failure hazards for each
observation at each time point in \code{valid_data}, if specified, or
\code{train_data} otherwise, had these observations been exposed. (3)
\code{"noexp_estimates"}, the expected conditional failure hazards for each
observation at each time point in \code{valid_data}, if specified, or
\code{train_data} otherwise, had these observations not been exposed.
}
\description{
\code{fit_failure_hazard()} estimates the conditional failure
hazard nuisance parameter. The estimator used for this estimation is based
on the \code{learners} argument, and the covariates considered are
specified by the \code{confounders} and \code{exposure} arguments.
}
\keyword{internal}
