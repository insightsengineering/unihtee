% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eifs.R
\name{uncentered_eif}
\alias{uncentered_eif}
\title{Uncentered Efficient Influence Function Computer}
\usage{
uncentered_eif(
  data,
  effect,
  confounders,
  exposure,
  outcome,
  modifiers,
  prop_score_fit,
  prop_score_values,
  cond_outcome_fit,
  failure_hazard_fit,
  censoring_hazard_fit
)
}
\arguments{
\item{data}{A \code{data.table} containing the observed data.
\code{train_data} is formatted by \code{\link{unihtee}()}.}

\item{effect}{A \code{character} indicating the type of treatment effect
modifier variable importance parameter. Currently supports
\code{"absolute"} and \code{"relative"}.}

\item{confounders}{A \code{character} vector of column names corresponding to
baseline covariates.}

\item{exposure}{A \code{character} corresponding to the exposure variable.}

\item{outcome}{A \code{character} corresponding to the outcome variable.}

\item{modifiers}{A \code{character} vector of columns names corresponding to
the suspected effect modifiers. This vector must be a subset of
\code{confounders}.}

\item{prop_score_fit}{A \code{list} output by the
\code{\link{fit_prop_score}()} function.}

\item{prop_score_values}{A \code{numeric} vector corresponding to the (known)
propensity score values for each observation in \code{data}.}

\item{cond_outcome_fit}{A \code{list} output by the
\code{\link{fit_failure_hazard}()} function.'}

\item{failure_hazard_fit}{A \code{list} output by the
\code{\link{fit_cond_outcome}()} function.'}

\item{censoring_hazard_fit}{A \code{list} output by the
\code{\link{fit_censoring_hazard}()} function.'}
}
\value{
A \code{data.table} whose columns are the uncentered efficient
influence functions of each variable in \code{modifiers}. The rows
correspond to the observations of \code{data}.
}
\description{
\code{uncentered_eif()} computes the efficient influence
function for the chosen parameter using the already estimated nuisance
parameters. If certain nuisance parameters are known, such as propensity
scores in a randomized control trial, then they may be input directly into
this function.
}
\keyword{internal}
