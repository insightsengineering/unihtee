utils::globalVariables(c("..to_keep", ".SD", "p_value"))
#' @title Univariate Heterogeneous Treatment Effect Modifier Estimator
#'
#' @description
#'
#' @details
#'
#' @param data
#' @param confounders
#' @param modifiers
#' @param exposure
#' @param outcome
#' @param cond_outcome_estimator
#' @param prop_score_estimator
#' @param prop_score_values
#'
#' @return
#'
#' @importFrom data.table as.data.table ".."
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
