% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unihtee.R
\name{cross_fit_fold}
\alias{cross_fit_fold}
\title{Cross-Fitting Procedure}
\usage{
cross_fit_fold(
  fold,
  data,
  confounders,
  modifiers,
  exposure,
  outcome,
  censoring,
  outcome_type,
  effect,
  estimator,
  cond_outcome_estimator,
  prop_score_estimator,
  prop_score_values,
  failure_hazard_estimator,
  censoring_hazard_estimator
)
}
\arguments{
\item{fold}{An \code{\link[origami:make_folds]{origami}} fold object.}

\item{data}{A \code{data.table} containing the observed data.}

\item{confounders}{A \code{character} vector of column names corresponding to
baseline covariates.}

\item{modifiers}{A \code{character} vector of columns names corresponding to
the suspected effect modifiers. This vector must be a subset of
\code{confounders}.}

\item{exposure}{A \code{character} corresponding to the exposure variable.}

\item{outcome}{A \code{character} corresponding to the outcome variable.}

\item{censoring}{A \code{character} indicating the right censoring indicator
variable. Only used with time-to-event outcomes. Defaults to \code{NULL}.}

\item{outcome_type}{A \code{character} indicating the outcome type.
\code{"continuous"}, \code{"binary"} and \code{"time-to-event"} are
currently supported.}

\item{effect}{A \code{character} indicating the type of treatment effect
modifier variable importance parameter. Currently supports
\code{"absolute"} and \code{"relative"}.}

\item{estimator}{A \code{character} set to either \code{"tmle"} or
\code{"onestep"}. The former results in \code{unihtee()} to use a targeted
maximum likelihood estimators to estimate the desired TEM-VIP, while the
latter uses a one step estimator.}

\item{cond_outcome_estimator}{A \code{\link[sl3]{Stack}}, or other learner
class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
learners from \pkg{sl3} to estimate the conditional outcome. Defaults to a
generalized linear model with one- and two- way interactions among all
\code{confounders} and \code{exposure} variables. Only used with continuous
and binary outcomes.}

\item{prop_score_estimator}{A \code{\link[sl3]{Stack}}, or other learner
class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
learners from \pkg{sl3} to estimate the propensity score. Defaults to a
generalized linear model with one- and two- way interactions among all
\code{confounders} variables.}

\item{prop_score_values}{An optional \code{character} corresponding to the
(known) propensity score values for each observation in \code{data}.
Defaults to \code{NULL}.}

\item{failure_hazard_estimator}{A \code{\link[sl3]{Stack}}, or other learner
class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
learners from \pkg{sl3} to estimate the conditional failure hazard
function. Defaults to an XGBoost learner with \code{confounders} and
\code{exposure} variables as covariates. Only used with time-to-event
outcomes.}

\item{censoring_hazard_estimator}{A \code{\link[sl3]{Stack}}, or other
learner class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a
set of learners from \pkg{sl3} to estimate the conditional censoring hazard
function. Defaults to an XGBoost learner with \code{confounders} and
\code{exposure} variables as covariates. Only used with time-to-event
outcomes.}
}
\value{
A \code{list} object containing the validation dataset's uncentered
efficient influence function estimates, the treatment effect modification
variable importance parameter estimates, the average causal effect
estimate, and the proportion of observations in the validation data.
}
\description{
\code{cross_fit_fold()} estimates the treatment effect
modification variable importance parameters using cross-fitting.
}
\keyword{internal}
