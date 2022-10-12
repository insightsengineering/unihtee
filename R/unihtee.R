utils::globalVariables(c("..to_keep", ".SD", "p_value"))
#' @title Univariate Heterogeneous Treatment Effect Modifier Estimator
#'
#' @description \code{unihtee()} estimates treatment effect modifiers variable
#'   importance parameters (TEM VIPs). These TEM VIPs are defined on the risk
#'   difference or relative risk scales and can be estimated using one step or
#'   targeted maximum likelihood estimators.
#'
#' @param data A \code{data.table} containing the observed data.
#' @param confounders A \code{character} vector of column names corresponding to
#'   baseline covariates.
#' @param modifiers A \code{character} vector of columns names corresponding to
#'   the suspected effect modifiers. This vector must be a subset of
#'   \code{confounders}.
#' @param exposure A \code{character} corresponding to the exposure variable.
#' @param outcome A \code{character} corresponding to the outcome variable.
#' @param cond_outcome_estimator A \code{\link[sl3]{Stack}}, or other learner
#'   class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
#'   learners from \pkg{sl3} to estimate the propensity score model. Defaults to
#'   a generalized linear model with one- and two- way interactions among all
#'   \code{confounders} and \code{exposure} variables.
#' @param prop_score_estimator A \code{\link[sl3]{Stack}}, or other learner
#'   class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
#'   learners from \pkg{sl3} to estimate the propensity score model. Defaults to
#'   a generalized linear model with one- and two- way interactions among all
#'   \code{confounders} variables.
#' @param prop_score_values A \code{numeric} vector corresponding to the (known)
#'   propensity score values for each observation in \code{data}.
#'
#' @return A \code{data.table} containing the effect estimates and (adjusted)
#'   p-values of the \code{modifiers}. The suspected treatment effect modifiers
#'   ordered according to ascending p-values.
#'
#' @importFrom data.table as.data.table

unihtee <- function(data,
                    confounders,
                    modifiers,
                    exposure,
                    outcome,
                    cond_outcome_estimator = sl3::Lrnr_glm_fast$new(),
                    prop_score_estimator = sl3::Lrnr_glm_fast$new(),
                    prop_score_values = NULL
                    ) {

  # transform data into data.table object with required columns
  data <- data.table::as.data.table(data)
  to_keep <- unique(c(confounders, modifiers, exposure, outcome))
  data <- data[, ..to_keep]

  # estimate nuisance parameters
  prop_score_fit <- fit_prop_score(
    train_data = data,
    valid_data = NULL,
    learners = prop_score_estimator,
    exposure = exposure,
    confounders = confounders
  )

  # fit the expected cond outcome
  cond_outcome_fit <- fit_cond_outcome(
    train_data = data,
    valid_data = NULL,
    learners = cond_outcome_estimator,
    exposure = exposure,
    confounders = confounders,
    outcome = outcome
  )

  # compute the efficient influence function
  ueif_dt <- uncentered_eif(
    data = data,
    confounders = confounders,
    exposure = exposure,
    outcome = outcome,
    modifiers = modifiers,
    prop_score_fit = prop_score_fit,
    prop_score_values = prop_score_values,
    cond_outcome_fit = cond_outcome_fit
  )

  # estimate the estimands
  one_step_fit <- one_step_estimator(uncentered_eif_data = ueif_dt)

  # compute the confidence intervals
  eif_vars <- ueif_dt[, lapply(.SD, var)]
  test_dt <- test_hypotheses(
    n_obs = nrow(data),
    estimates = one_step_fit,
    var_estimates = eif_vars
  )

  # organize table in decreasing order of p value
  test_dt <- test_dt[order(p_value), ]

  return(test_dt)

}
