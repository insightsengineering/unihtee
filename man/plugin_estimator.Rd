% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimators.R
\name{plugin_estimator}
\alias{plugin_estimator}
\title{Plug-In Estimator}
\usage{
plugin_estimator(
  data,
  outcome,
  modifiers,
  effect,
  cond_outcome_fit,
  failure_hazard_fit
)
}
\arguments{
\item{data}{A \code{data.table} containing the observed data.
\code{train_data} is formatted by \code{\link{unihtee}()}.}

\item{outcome}{A \code{character} corresponding to the outcome variable.}

\item{modifiers}{A \code{character} vector of columns names corresponding to
the suspected effect modifiers. This vector must be a subset of
\code{confounders}.}

\item{effect}{A \code{character} indicating the type of treatment effect
modifier variable importance parameter. Currently supports
\code{"absolute"} and \code{"relative"}.}

\item{cond_outcome_fit}{A \code{list} output by the
\code{\link{fit_cond_outcome}()} function.}

\item{failure_hazard_fit}{A \code{list} output by the
\code{\link{fit_cond_outcome}()} function.'}
}
\value{
A vector containing the plug-in estimates for each potential
modifier.
}
\description{
\code{plugin_estimator()} estimates the treatment effect modifier variable
importance parameters using a plug-in estimation strategy.
}
\keyword{internal}
