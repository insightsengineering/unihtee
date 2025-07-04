#' @title Propensity Score Estimator
#'
#' @description \code{fit_prop_score()} estimates the propensity score over the
#'   \code{train_data} dataset. The estimator used for this estimation is based
#'   on the \code{learners} argument, and the covariates considered are
#'   specified by the \code{confounders} argument.
#'
#' @param train_data A \code{data.table} containing the observed data.
#'   \code{train_data} is formatted by \code{\link{unihtee}()}.
#' @param valid_data An optional \code{data.table} representing a holdout
#'   dataset of the observed data. It is only used for cross-fitting purposes.
#'   Defaults to \code{NULL}.
#' @param learners A \code{\link[sl3]{Stack}}, or other learner class
#'   (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
#'   learners from \pkg{sl3} to estimate the propensity score model.
#' @param exposure  A \code{character} corresponding to the exposure variable.
#' @param confounders A \code{character} vector of column names corresponding to
#'   baseline covariates.
#'
#' @return A named \code{list} of two elements. (1) \code{"estimates"}, the
#'   propensity score estimates for each observation in \code{valid_data}, if
#'   specified, or \code{train_data} otherwise. (2) \code{"fit"}, the trained
#'   \code{\link[sl3]{Stack}} or learner.
#'
#' @importFrom sl3 sl3_Task
#'
#' @keywords internal
fit_prop_score <- function(train_data,
                           valid_data,
                           learners,
                           exposure,
                           confounders) {

  # construct the training task
  train_data_task <- sl3::sl3_Task$new(
    data = train_data,
    covariates = confounders,
    outcome = exposure,
    outcome_type = "binomial"
  )

  # estimate the propensity score
  prop_score_fit <- learners$train(train_data_task)

  if (is.null(valid_data)) {

    # extract the propensity score predictions for the training data
    prop_score_est <- prop_score_fit$predict()
  } else {

    # construct trask for validation dataset prediction
    valid_data_task <- sl3::sl3_Task$new(
      data = valid_data,
      covariates = confounders,
      outcome = exposure,
      outcome_type = "binomial"
    )

    # compute the propensity score estimates on the validation dataset
    prop_score_est <- prop_score_fit$predict(valid_data_task)
  }


  ## bound the propensity score estimates to avoid practical positivity
  ## violations
  truncation <- 1 / min(sqrt(nrow(train_data)), 1000)
  prop_score_est[prop_score_est < truncation] <- truncation
  prop_score_est[prop_score_est > (1 - truncation)] <- (1 - truncation)

  return(list(
    "estimates" = prop_score_est,
    "fit" = prop_score_fit
  ))
}

#' @title Conditional Outcome Estimator
#'
#' @description \code{fit_cond_outcome()} estimates the conditional outcome,
#'   based on the exposure and confounders, over the \code{train_data} dataset.
#'   The estimator used for this estimation is based on the \code{learners}
#'   argument, and the covariates considered are specified by the
#'   \code{confounders} and \code{exposure} arguments.
#'
#' @param train_data A \code{data.table} containing the observed data.
#'   \code{train_data} is formatted by \code{\link{unihtee}()}.
#' @param valid_data An optional \code{data.table} representing a holdout
#'   dataset of the observed data. It is only used for cross-fitting purposes.
#'   Defaults to \code{NULL}.
#' @param learners A \code{\link[sl3]{Stack}}, or other learner class
#'   (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
#'   learners from \pkg{sl3} to estimate the propensity score model.
#' @param outcome  A \code{character} corresponding to the outcome variable.
#' @param exposure  A \code{character} corresponding to the exposure variable.
#' @param confounders A \code{character} vector of column names corresponding to
#'   baseline covariates.
#'
#' @return A named \code{list} of four elements. (1) \code{"estimates"}, the
#'   expected conditional outcome for each observation in \code{valid_data}, if
#'   specified, or \code{train_data} otherwise. (2) \code{"fit"}, the trained
#'   \code{\link[sl3]{Stack}} or learner. (3) \code{"exp_estimates"}, the
#'   expected conditional outcome for each observation in \code{valid_data}, if
#'   specified, or \code{train_data} otherwise, had these observations been
#'   exposed. (4) \code{"noexp_estimates"}, the expected conditional outcome for
#'   each observation in \code{valid_data}, if specified, or \code{train_data}
#'   otherwise, had these observations not been exposed.
#'
#' @importFrom sl3 sl3_Task
#' @importFrom data.table copy
#'
#' @keywords internal
fit_cond_outcome <- function(train_data,
                             valid_data,
                             learners,
                             outcome,
                             exposure,
                             confounders) {

  # define the covariates
  covariates <- c(exposure, confounders)

  # construct the training task
  train_data_task <- sl3::sl3_Task$new(
    data = train_data,
    covariates = covariates,
    outcome = outcome
  )

  # estimate the conditional outcome
  cond_outcome_fit <- learners$train(train_data_task)

  if (is.null(valid_data)) {

    # extract the cond outcome predictions for the training data
    cond_outcome_est <- cond_outcome_fit$predict()

    # copy the train data to the potential outcome datasets
    exp_data <- data.table::copy(train_data)
    noexp_data <- data.table::copy(train_data)
  } else {

    # construct the validation data task
    valid_data_task <- sl3::sl3_Task$new(
      data = valid_data,
      covariates = covariates,
      outcome = outcome
    )

    # extract the cond outcome predictions for the valid data
    cond_outcome_est <- cond_outcome_fit$predict(valid_data_task)

    # copy the train data to the potential outcome datasets
    exp_data <- data.table::copy(valid_data)
    noexp_data <- data.table::copy(valid_data)
  }

  # estimate the potential outcomes
  exp_data[[exposure]] <- 1
  exp_data_task <- sl3::sl3_Task$new(
    data = exp_data,
    covariates = covariates,
    outcome = outcome
  )
  exp_cond_outcome_est <- cond_outcome_fit$predict(exp_data_task)
  noexp_data[[exposure]] <- 0
  noexp_data_task <- sl3::sl3_Task$new(
    data = noexp_data,
    covariates = covariates,
    outcome = outcome
  )
  noexp_cond_outcome_est <- cond_outcome_fit$predict(noexp_data_task)


  return(list(
    "estimates" = cond_outcome_est,
    "fit" = cond_outcome_fit,
    "exp_estimates" = exp_cond_outcome_est,
    "noexp_estimates" = noexp_cond_outcome_est
  ))
}

#' @title Conditional Failure Hazard Estimator
#'
#' @description \code{fit_failure_hazard()} estimates the conditional failure
#'   hazard nuisance parameter. The estimator used for this estimation is based
#'   on the \code{learners} argument, and the covariates considered are
#'   specified by the \code{confounders} and \code{exposure} arguments.
#'
#' @param train_data A long \code{data.table} containing the observed data.
#'   \code{train_data} is formatted by \code{\link{unihtee}()} and
#'   \code{\link{tte_data_melt}()}.
#' @param valid_data An optional \code{data.table} representing a holdout
#'   dataset of the observed data. It is only used for cross-fitting purposes.
#'   Defaults to \code{NULL}.
#' @param learners A \code{\link[sl3]{Stack}}, or other learner class
#'   (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
#'   learners from \pkg{sl3} to estimate the propensity score model.
#' @param exposure A \code{character} corresponding to the exposure variable.
#' @param confounders A \code{character} vector of column names corresponding to
#'   baseline covariates.
#' @param times A \code{character} giving the name corresponding to the times.
#'
#' @return A named \code{list} of three elements. (1) \code{"estimates"}, the
#'   expected failure hazards for each observation at each time point in
#'   \code{valid_data}, if specified, or \code{train_data} otherwise. (2)
#'   \code{"exp_estimates"}, the expected conditional failure hazards for each
#'   observation at each time point in \code{valid_data}, if specified, or
#'   \code{train_data} otherwise, had these observations been exposed. (3)
#'   \code{"noexp_estimates"}, the expected conditional failure hazards for each
#'   observation at each time point in \code{valid_data}, if specified, or
#'   \code{train_data} otherwise, had these observations not been exposed.
#'
#' @importFrom sl3 sl3_Task
#' @importFrom data.table copy
#'
#' @keywords internal
#'
fit_failure_hazard <- function(train_data,
                               valid_data,
                               learners,
                               confounders,
                               exposure,
                               times) {

  # define the covariates
  covariates <- c(times, exposure, confounders)

  # construct the training task
  keep_var <- "keep"
  failure_hazard_task <- sl3::sl3_Task$new(
    data = train_data[get(keep_var) == 1, ],
    outcome = "failure",
    covariates = covariates
  )

  # estimate the conditional failure hazard function
  failure_hazard_fit <- learners$train(failure_hazard_task)

  # compute estimates
  if (is.null(valid_data)) {

    # extract the cond failure hazard predictions for the training data
    train_failure_hazard_task <- sl3::sl3_Task$new(
      data = train_data,
      outcome = "failure",
      covariates = covariates
    )
    estimates <- failure_hazard_fit$predict(train_failure_hazard_task)

    # copy the train data to the potential outcome datasets
    exp_data <- data.table::copy(train_data)
    noexp_data <- data.table::copy(train_data)

  } else {

    # construct the validation data task
    valid_data_task <- sl3::sl3_Task$new(
      data = valid_data,
      covariates = covariates,
      outcome = "failure"
    )

    # extract the cond failure hazard predictions for the valid data
    estimates <- failure_hazard_fit$predict(valid_data_task)

    # copy the train data to the potential outcome datasets
    exp_data <- data.table::copy(valid_data)
    noexp_data <- data.table::copy(valid_data)

  }

  # estimate the failure hazards under each exposure level
  exp_data[[exposure]] <- 1
  exp_data_task <- sl3::sl3_Task$new(
    data = exp_data,
    covariates = covariates,
    outcome = "failure"
  )
  exp_estimates <- failure_hazard_fit$predict(exp_data_task)
  noexp_data[[exposure]] <- 0
  noexp_data_task <- sl3::sl3_Task$new(
    data = noexp_data,
    covariates = covariates,
    outcome = "failure"
  )
  noexp_estimates <- failure_hazard_fit$predict(noexp_data_task)


  return(list(
    "estimates" = estimates,
    "exp_estimates" = exp_estimates,
    "noexp_estimates" = noexp_estimates
  ))
}

#' @title Conditional Censoring Hazard Estimator
#'
#' @description \code{fit_censoring_hazard()} estimates the conditional
#'   censoring hazard nuisance parameter. The estimator used for this estimation
#'   is based on the \code{learners} argument, and the covariates considered are
#'   specified by the \code{confounders} and \code{exposure} arguments.
#'
#' @param train_data A long \code{data.table} containing the observed data.
#'   \code{train_data} is formatted by \code{\link{unihtee}()} and
#'   \code{\link{tte_data_melt}()}.
#' @param valid_data An optional \code{data.table} representing a holdout
#'   dataset of the observed data. It is only used for cross-fitting purposes.
#'   Defaults to \code{NULL}.
#' @param learners A \code{\link[sl3]{Stack}}, or other learner class
#'   (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
#'   learners from \pkg{sl3} to estimate the propensity score model.
#' @param exposure A \code{character} corresponding to the exposure variable.
#' @param confounders A \code{character} vector of column names corresponding to
#'   baseline covariates.
#' @param times A \code{character} giving the name corresponding to the times.
#' @param censoring A \code{character} indicating the censoring indicator.
#'
#' @return A named \code{list} of three elements. (1) \code{"estimates"}, the
#'   expected censoring hazards for each observation at each time point in
#'   \code{valid_data}, if specified, or \code{train_data} otherwise. (2)
#'   \code{"exp_estimates"}, the expected conditional censoring hazards for each
#'   observation at each time point in \code{valid_data}, if specified, or
#'   \code{train_data} otherwise, had these observations been exposed. (3)
#'   \code{"noexp_estimates"}, the expected conditional censoring hazards for
#'   each observation at each time point in \code{valid_data}, if specified, or
#'   \code{train_data} otherwise, had these observations not been exposed.
#'
#' @importFrom sl3 sl3_Task
#' @importFrom data.table copy
#'
#' @keywords internal
#'
fit_censoring_hazard <- function(train_data,
                                 valid_data,
                                 learners,
                                 confounders,
                                 exposure,
                                 times,
                                 censoring) {

  # define the covariates
  covariates <- c(times, exposure, confounders)

  # construct the training task
  keep_var <- "keep"
  censoring_hazard_task <- sl3::sl3_Task$new(
    data = train_data[get(keep_var) == 1],
    outcome = censoring,
    covariates = covariates
  )

  # estimate the conditional censoring hazard function
  censoring_hazard_fit <- learners$train(censoring_hazard_task)

  # compute estimates
  if (is.null(valid_data)) {

    # extract the cond censoring hazard predictions for the training data
    full_censoring_hazard_task <- sl3::sl3_Task$new(
      data = train_data,
      outcome = censoring,
      covariates = covariates
    )
    estimates <- censoring_hazard_fit$predict(full_censoring_hazard_task)

    # copy the train data to the potential outcome datasets
    exp_data <- data.table::copy(train_data)
    noexp_data <- data.table::copy(train_data)

  } else {

    # construct the validation data task
    valid_data_task <- sl3::sl3_Task$new(
      data = valid_data,
      covariates = covariates,
      outcome = censoring
    )

    # extract the cond censoring hazard predictions for the valid data
    estimates <- censoring_hazard_fit$predict(valid_data_task)

    # copy the train data to the potential outcome datasets
    exp_data <- data.table::copy(valid_data)
    noexp_data <- data.table::copy(valid_data)

  }

  # estimate the censoring hazards under each exposure level
  exp_data[[exposure]] <- 1
  exp_data_task <- sl3::sl3_Task$new(
    data = exp_data,
    covariates = covariates,
    outcome = censoring
  )
  exp_estimates <- censoring_hazard_fit$predict(exp_data_task)
  noexp_data[[exposure]] <- 0
  noexp_data_task <- sl3::sl3_Task$new(
    data = noexp_data,
    covariates = covariates,
    outcome = censoring
  )
  noexp_estimates <- censoring_hazard_fit$predict(noexp_data_task)


  return(list(
    "estimates" = estimates,
    "exp_estimates" = exp_estimates,
    "noexp_estimates" = noexp_estimates
  ))
}


one_step_ate_estimator <- function(
  data,
  confounders,
  exposure,
  outcome,
  prop_score_fit,
  prop_score_values = NULL,
  cond_outcome_fit
) {

  ## compute the inverse probability weights
  if (!is.null(prop_score_values)) {
    prop_scores <- data[[prop_score_values]]
  } else {
    prop_scores <- prop_score_fit$estimates
  }
  ipws <- (2 * data[[exposure]] - 1) /
    (data[[exposure]] * prop_scores +
       (1 - data[[exposure]]) * (1 - prop_scores))

  ## compute conditional outcome residuals
  cond_outcome_resid <- data[[outcome]] - cond_outcome_fit$estimates

  # compute the uncentered EIF of the ATE
  uncentered_eif <- ipws * cond_outcome_resid + cond_outcome_fit$exp_estimates -
    cond_outcome_fit$noexp_estimates

  # compute the estimate
  estimate <- mean(uncentered_eif)

  return(estimate)
}

tml_ate_estimator <- function(
    data,
    confounders,
    exposure,
    outcome,
    prop_score_fit,
    prop_score_values = NULL,
    cond_outcome_fit
) {

  ## compute the inverse probability weights
  if (!is.null(prop_score_values)) {
    prop_scores <- data[[prop_score_values]]
  } else {
    prop_scores <- prop_score_fit$estimates
  }

  ## tilt conditional expected outcome under exposure
  exp_estimates <- cond_outcome_fit$exp_estimates
  h_num_1 <- data[[exposure]]
  h_denom_1 <- data[[exposure]] * prop_scores
  q_1_star_logit <- stats::qlogis(bound_precision(exp_estimates))
  suppressWarnings(
    q_1_tilt_fit <- stats::glm(
      data[[outcome]] ~ -1 + h_num_1,
      offset = q_1_star_logit,
      weights = h_denom_1,
      family = "binomial"
    )
  )
  q_1_star <- predict(q_1_tilt_fit, type = "response")

  ## tilt conditional expected outcome under non-exposure
  noexp_estimates <- cond_outcome_fit$noexp_estimates
  h_num_0 <- 1 - data[[exposure]]
  h_denom_0 <- (1 - data[[exposure]]) / (1 - prop_scores)
  q_0_star_logit <- stats::qlogis(bound_precision(noexp_estimates))
  suppressWarnings(
    q_0_tilt_fit <- stats::glm(
      data[[outcome]] ~ -1 + h_num_0,
      offset = q_0_star_logit,
      weights = h_denom_0,
      family = "binomial"
    )
  )
  q_0_star <- predict(q_0_tilt_fit, type = "response")

  ## compute the plug-in estimate
  estimate <- mean(q_1_star - q_0_star)

  return(estimate)

}
