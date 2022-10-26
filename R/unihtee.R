utils::globalVariables(c("..to_keep", ".SD", "p_value"))
#' @title Univariate Heterogeneous Treatment Effect Modifier Estimator
#'
#' @description \code{unihtee()} estimates treatment effect modifiers variable
#'   importance parameters (TEM VIPs). These TEM VIPs are defined on the risk
#'   difference or relative risk scales and can be estimated using one-step or
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
#' @param censoring A \code{character} indicating the right censoring indicator
#'   variable. Only used with time-to-event outcomes. Defaults to \code{NULL}.
#' @param time_cutoff A \code{numeric} representing the timepoint at which to
#'   evaluate the time-to-event parameter. Only used with time-to-event
#'   outcomes. Defaults to \code{NULL}.
#' @param outcome_type A \code{character} indicating the outcome type.
#'   \code{"continuous"}, \code{"binary"} and \code{"time-to-event"} are
#'   currently supported.
#' @param risk_type A \code{character} indicating the type of treatment effect
#'   modifier variable importance parameter. Currently supports
#'   \code{"risk difference"} and \code{"relative risk"}.
#' @param estimator A \code{character} set to either \code{"tmle"} or
#'   \code{"onestep"}. The former results in \code{unihtee()} to use a targeted
#'   maximum likelihood estimators to estimate the deisred TEM VIP, while the
#'   latter uses a one step estimator.
#' @param cond_outcome_estimator A \code{\link[sl3]{Stack}}, or other learner
#'   class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
#'   learners from \pkg{sl3} to estimate the conditional outcome. Defaults to a
#'   generalized linear model with one- and two- way interactions among all
#'   \code{confounders} and \code{exposure} variables. Only used with continuous
#'   and binary outcomes.
#' @param prop_score_estimator A \code{\link[sl3]{Stack}}, or other learner
#'   class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
#'   learners from \pkg{sl3} to estimate the propensity score. Defaults to a
#'   generalized linear model with one- and two- way interactions among all
#'   \code{confounders} variables.
#' @param prop_score_values An optional \code{character} corresponding to the
#'   (known) propensity score values for each observation in \code{data}.
#'   Defaults to \code{NULL}.
#' @param failure_hazard_estimator A \code{\link[sl3]{Stack}}, or other learner
#'   class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
#'   learners from \pkg{sl3} to estimate the conditional failure hazard
#'   function. Defaults to an xgboost learner with \code{confounders} and
#'   \code{exposure} variables as covariates. Only used with time-to-event
#'   outcomes.
#' @param censoring_hazard_estimator A \code{\link[sl3]{Stack}}, or other
#'   learner class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a
#'   set of learners from \pkg{sl3} to estimate the conditional censoring hazard
#'   function. Defaults to an xgboost learner with \code{confounders} and
#'   \code{exposure} variables as covariates. Only used with time-to-event
#'   outcomes.
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
                    censoring = NULL,
                    time_cutoff = NULL,
                    outcome_type = c("continuous", "binary", "time-to-event"),
                    risk_type = c("risk difference", "relative risk"),
                    estimator = c("tmle", "onestep"),
                    cond_outcome_estimator = sl3::Lrnr_glm_fast$new(),
                    prop_score_estimator = sl3::Lrnr_glm_fast$new(),
                    prop_score_values = NULL,
                    failure_hazard_estimator = sl3::Lrnr_xgboost$new(),
                    censoring_hazard_estimator = sl3::Lrnr_xgboost$new()) {

  ## specify the TEM VIP type
  risk_type <- match.arg(risk_type)

  ## specify the outcome type
  outcome_type <- match.arg(outcome_type)

  ## set the estimator
  estimator <- match.arg(estimator)

  ## transform data into data.table object with required columns
  data <- data.table::as.data.table(data)
  to_keep <- unique(c(
    confounders, modifiers, exposure, outcome, prop_score_values, censoring
  ))
  data <- data[, ..to_keep]

  ## scale the outcome to be between 0 and 1 if outcome is continuous
  if (outcome_type == "continuous") {
    min_out <- min(data[[outcome]])
    max_out <- max(data[[outcome]])
    rescale_factor <- max_out - min_out
    data[[outcome]] <- (data[[outcome]] - min_out) / rescale_factor
  } else {
    rescale_factor <- 1
  }

  ## create the long-form data for time-to-event outcomes
  if (outcome_type == "time-to-event") {
    long_data <- tte_data_melt(
      data = data,
      confounders = confounders,
      exposure = exposure,
      outcome = outcome,
      censoring = censoring,
      time_cutoff = time_cutoff,
      prop_score_values = prop_score_values
    )
  }

  ## estimate nuisance parameters

  ## fit the propensity score, if necessary
  if (is.null(prop_score_values)) {
    prop_score_fit <- fit_prop_score(
      train_data = data,
      valid_data = if (outcome_type == "time-to-event") long_data else NULL,
      learners = prop_score_estimator,
      exposure = exposure,
      confounders = confounders
    )
  } else {
    prop_score_fit <- list(
      estimates = data[[prop_score_values]]
    )
  }

  ## fit the remaining nuisance parameter for continuous and binary outcomes
  if (outcome_type %in% c("continuous", "binary")) {
    cond_outcome_fit <- fit_cond_outcome(
      train_data = data,
      valid_data = NULL,
      learners = cond_outcome_estimator,
      exposure = exposure,
      confounders = confounders,
      outcome = outcome
    )
    failure_hazard_fit <- NULL
    censoring_hazard_fit <- NULL
  }

  ## fit the remaining nuisance paramters for time-to-event outcomes
  if (outcome_type == "time-to-event") {
    failure_hazard_fit <- fit_failure_hazard(
      train_data = long_data,
      valid_data = NULL,
      learners = failure_hazard_estimator,
      exposure = exposure,
      confounders = confounders
    )
    censoring_hazard_fit <- fit_censoring_hazard(
      train_data = long_data,
      valid_data = NULL,
      learners = censoring_hazard_estimator,
      exposure = exposure,
      confounders = confounders,
      censoring = censoring
    )
    cond_outcome_fit <- NULL
  }

  ## compute the efficient influence function
  ueif_dt <- uncentered_eif(
    data = if (outcome_type == "time-to-event") long_data else data,
    type = risk_type,
    confounders = confounders,
    exposure = exposure,
    outcome = outcome,
    modifiers = modifiers,
    prop_score_fit = prop_score_fit,
    prop_score_values = prop_score_values,
    cond_outcome_fit = cond_outcome_fit,
    failure_hazard_fit = failure_hazard_fit,
    censoring_hazard_fit = censoring_hazard_fit
  )

  ## estimate the estimands
  if (estimator == "onestep") {
    tem_vip_fit <- one_step_estimator(uncentered_eif_data = ueif_dt)
  } else {
    tem_vip_fit <- tml_estimator(
      data = if (outcome_type == "time-to-event") long_data else data,
      confounders = confounders,
      modifiers = modifiers,
      exposure = exposure,
      outcome = outcome,
      type = risk_type,
      prop_score_fit = prop_score_fit,
      cond_outcome_fit = cond_outcome_fit,
      failure_hazard_fit = failure_hazard_fit,
      censoring_hazard_fit = censoring_hazard_fit
    )
  }

  ## compute the confidence intervals and rescale everything
  eif_vars <- ueif_dt[, lapply(.SD, var)]
  test_dt <- test_hypotheses(
    n_obs = nrow(data),
    estimates = tem_vip_fit,
    var_estimates = eif_vars,
    rescale_factor = rescale_factor
  )

  ## organize table in decreasing order of p-value
  test_dt <- test_dt[order(p_value), ]

  return(test_dt)
}
