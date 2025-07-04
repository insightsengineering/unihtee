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


#' One Step Estimator of the Average Treatment Effect
#'
#' `one_step_ate_estimator()` implements a one-step estimator of the average
#' treatment effect.
#'
#' @inheritParams tml_estimator
#'
#' @returns The one-step estimate of the average treatment effect.
#'
#' @keywords internal
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

#' Targeted Maximum Likelihood Estimator of the Average Treatment Effect
#'
#' `tml_ate_estimator()` implements a targeted maximum likelihood estimator of
#' the average treatment effect.
#'
#' @inheritParams tml_estimator
#'
#' @returns The targeted maximum likelihood estimate of the average treatment
#'   effect.
#'
#' @keywords internal
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

#' One-Step Estimator of the Marginal Relative Effect for Non-Time-to-Event Outcomes
#'
#' @inheritParams tml_estimator
#'
#' @returns The one-step estimate.
#'
#' @keywords internal
one_step_estimator_ate_log_outcome <- function(
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

  # compute the uncentered EIF
  eps <- 1e-10
  estimates <- cond_outcome_fit$estimates
  estimates[estimates < eps] <- eps
  exp_estimates <- cond_outcome_fit$exp_estimates
  exp_estimates[exp_estimates < eps] <- eps
  noexp_estimates <- cond_outcome_fit$noexp_estimates
  noexp_estimates[noexp_estimates < eps] <- eps
  aipws <- ipws * cond_outcome_resid / estimates +
    log(exp_estimates) - log(noexp_estimates)

  # compute the estimate
  estimate <- mean(aipws)

  return(estimate)
}

#' TMLE of the Marginal Relative Effect for Non-Time-to-Event Outcomes
#'
#' @inheritParams tml_estimator
#'
#' @returns The TML estimate.
#'
#' @keywords internal
tml_estimator_ate_log_outcome <- function(
    data,
    confounders,
    exposure,
    outcome,
    prop_score_fit,
    prop_score_values = NULL,
    cond_outcome_fit
) {

  # constant used to truncate estimates near zero
  eps <- .Machine$double.eps

  ## compute the inverse probability weights
  if (!is.null(prop_score_values)) {
    prop_scores <- data[[prop_score_values]]
  } else {
    prop_scores <- prop_score_fit$estimates
  }
  ipws <- (2 * data[[exposure]] - 1) /
    (data[[exposure]] * prop_scores +
       (1 - data[[exposure]]) * (1 - prop_scores))

  # truncate cond outcome estimates if necessary
  estimates <- cond_outcome_fit$estimates
  exp_estimates <- cond_outcome_fit$exp_estimates
  noexp_estimates <- cond_outcome_fit$noexp_estimates
  exp_estimates[exp_estimates < eps] <- eps
  noexp_estimates[noexp_estimates < eps] <- eps

  # compute that partial clever covariate
  h_num_1 <- data[[exposure]]
  h_denom_1 <- 1 / prop_scores
  h_num_0 <- (1 - data[[exposure]])
  h_denom_0 <- 1 / (1 - prop_scores)

  ## set the TML tilting hyperparameters
  score_stop_crit <- 1 / nrow(data) * log(nrow(data))
  tilt_tol <- 5
  max_iter <- 10

  ## compute the tilted conditional outcome estimators
  q_score <- Inf
  iter <- 1
  q_star <- estimates
  q_1_star <- exp_estimates
  q_0_star <- noexp_estimates

  while (score_stop_crit < abs(mean(q_score)) && iter < max_iter) {

    ## transform expected conditional outcomes for tilting
    q_star_logit <- stats::qlogis(bound_precision(q_star))
    q_1_star_logit <- stats::qlogis(bound_precision(q_1_star))
    q_0_star_logit <- stats::qlogis(bound_precision(q_0_star))

    ## tilt the clever covariate using the log loss for q_1
    suppressWarnings(
      q_1_tilt_fit <- stats::glm(
        data[[outcome]] ~ -1 + h_num_1,
        offset = q_1_star_logit,
        weights = h_denom_1,
        family = "binomial"
      )
    )

    ## tilt the clever covariate using the log loss for q_0
    suppressWarnings(
      q_0_tilt_fit <- stats::glm(
        data[[outcome]] ~ -1 + h_num_0,
        offset = q_0_star_logit,
        weights = h_denom_0,
        family = "binomial"
      )
    )

    ## update the nuisance parameter estimates
    q_1_star <- predict(q_1_tilt_fit, type = "response")
    q_0_star <- predict(q_0_tilt_fit, type = "response")
    q_star <- data[[exposure]] * q_1_star +
      (1 - data[[exposure]]) * q_0_star

    # avoid division by near zero in the score
    bounded_q_star <- q_star
    bounded_q_star[bounded_q_star < 0.1] <- 0.1

    ## compute the score
    q_score <- (2 * data[[exposure]] - 1) /
      (data[[exposure]] * prop_scores +
         (1 - data[[exposure]]) * (1 - prop_scores)) *
      (data[[outcome]] - q_star) / bounded_q_star

    iter <- iter + 1

  }

  ## compute the plugin estimate with the update cond outcome estimates
  estimate <- mean(log(q_1_star / q_0_star))

  return(estimate)
}


one_step_rmst_diff_estimator <- function(
  data,
  exposure,
  outcome,
  prop_score_fit,
  prop_score_values,
  failure_hazard_fit,
  censoring_hazard_fit
) {

  ## compute the inverse probability weights
  ## NOTE: PS estimates must be as long as long_dt in TTE settings
  if (!is.null(prop_score_values)) {
    prop_scores <- data[[prop_score_values]]
  } else {
    prop_scores <- prop_score_fit$estimates
  }
  ipws <- (2 * data[[exposure]] - 1) /
    (data[[exposure]] * prop_scores +
       (1 - data[[exposure]]) * (1 - prop_scores))

  ## compute the survival probabilities
  data$ipws <- -ipws
  data$failure_haz_est <- failure_hazard_fit$estimates
  data$failure_haz_exp_est <- failure_hazard_fit$exp_estimates
  data$failure_haz_noexp_est <- failure_hazard_fit$noexp_estimates
  data$censoring_haz_est <- censoring_hazard_fit$estimates
  data[, surv_est := cumprod(1 - failure_haz_est), by = "id"]
  data[, `:=`(
    surv_exp_est = cumprod(1 - failure_haz_exp_est),
    surv_noexp_est = cumprod(1 - failure_haz_noexp_est),
    cens_est = cumprod(1 - censoring_haz_est)
  ), by = "id"]
  data[, cens_est_lag := data.table::shift(cens_est, n = 1, fill = 1),
       by = "id"]

  ## compute the integrand of the EIF for each row
  data[, prev_time := data.table::shift(get(outcome), n = 1, fill = 0),
       by = "id"]
  data[, int_weight := as.numeric(get(outcome)) - as.numeric(prev_time),
       by = "id"]
  data[, inner_integrand := int_weight * keep * ipws /
         (cens_est_lag * surv_est) * (failure - failure_haz_est),
       by = "id"
  ]
  data[, inner_integral := cumsum(inner_integrand), by = "id"]
  data[, outer_integrand := int_weight * (surv_est * inner_integral +
                                            surv_exp_est - surv_noexp_est),
       by = "id"
  ]
  data[, aipws := cumsum(outer_integrand), by = "id"]

  ## use only the observations at the time_cutoff
  time_cutoff <- max(data[[outcome]])
  data <- data[get(outcome) == time_cutoff, ]

  ## return the AIPWs at the currect cutoff
  aipws <- data$aipws

  ## compute the one-step estimate
  estimate <- mean(aipws)

  return(estimate)
}

tml_rmst_diff_estimator <- function(
    data,
    exposure,
    outcome,
    prop_score_fit,
    prop_score_values,
    failure_hazard_fit,
    censoring_hazard_fit
) {

  # constant used to truncate estimates near zero
  eps <- .Machine$double.eps

  ## compute the inverse probability weights
  ## NOTE: PS estimates must be as long as long_dt in TTE settings
  if (!is.null(prop_score_values)) {
    prop_scores <- data[[prop_score_values]]
  } else {
    prop_scores <- prop_score_fit$estimates
  }

  ## compute the IPWs
  data$prop_scores <- prop_scores
  data$ipws <- -(2 * data[[exposure]] - 1) /
    (data[[exposure]] * prop_scores +
       (1 - data[[exposure]]) * (1 - prop_scores))

  ## compute the survival and censoring probabilities
  data$failure_haz_est <- failure_hazard_fit$estimates
  data$failure_haz_exp_est <- failure_hazard_fit$exp_estimates
  data$failure_haz_noexp_est <- failure_hazard_fit$noexp_estimates
  data$censoring_haz_est <- censoring_hazard_fit$estimates
  data[, surv_est := cumprod(1 - failure_haz_est), by = "id"]
  data[, `:=`(
    surv_exp_est = cumprod(1 - failure_haz_exp_est),
    surv_noexp_est = cumprod(1 - failure_haz_noexp_est),
    cens_est = cumprod(1 - censoring_haz_est)
  ), by = "id"]
  data[, cens_est_lag := data.table::shift(cens_est, n = 1, fill = 1),
       by = "id"]

  # compute the weights for envetual discrete integration
  data[, prev_time := data.table::shift(get(outcome), n = 1, fill = 0),
       by = "id"]
  data[, int_weight := as.numeric(get(outcome)) - as.numeric(prev_time),
       by = "id"]

  ## compute the partial clever covariate at each timepoint
  data[, inner_integrand := keep * ipws / (cens_est_lag * surv_est),
       by = "id"]

  ## set parameters for update step
  max_iter <- 10
  tilt_tol <- 5

  ## compute the partial clever covariates at each timepoint
  the_times <- sort(unique(data[[outcome]]))
  survival_preds <- lapply(
    the_times,
    function(t) {

      ## extract timepoints used for tilting
      filtered_dt <- data.table::copy(data)
      filtered_dt <- filtered_dt[get(outcome) <= t, ]

      ## the time cutoff
      filtered_dt[, surv_time_cutoff := min(surv_est)]

      ## compute the partial clever covariates
      filtered_dt[, `:=`(
        partial_h = inner_integrand * surv_time_cutoff,
        partial_h_1 = -surv_time_cutoff * keep /
          (cens_est_lag * surv_est),
        partial_h_0 = -surv_time_cutoff * keep /
          (cens_est_lag * surv_est)
      )]

      ## finalize the clever covariate
      mod_h <- filtered_dt$partial_h
      mod_h_1 <- filtered_dt$partial_h_1
      mod_h_0 <- filtered_dt$partial_h_0

      ## set initial values
      fail_haz_est <- filtered_dt$failure_haz_est
      fail_haz_1_est <- filtered_dt$failure_haz_exp_est
      fail_haz_0_est <- filtered_dt$failure_haz_noexp_est
      iter <- 1
      haz_score <- Inf
      score_stop_crit <- 1 /
        (nrow(filtered_dt) * log(nrow(filtered_dt)))

      ## tilt the conditional failure estimates
      while (abs(mean(haz_score)) > score_stop_crit &&
             iter < max_iter) {

        ## transform expected conditional outcomes for tilting
        fail_haz_est_logit <- stats::qlogis(
          bound_precision(fail_haz_est)
        )
        fail_haz_1_est_logit <- stats::qlogis(
          bound_precision(fail_haz_1_est)
        )
        fail_haz_0_est_logit <- stats::qlogis(
          bound_precision(fail_haz_0_est)
        )

        ## tilt the clever covariate of failure hazard under treatment
        suppressWarnings(
          fail_haz_1_tilt_fit <- stats::glm(
            filtered_dt$failure ~ -1 + mod_h_1,
            offset = fail_haz_1_est_logit,
            weights = filtered_dt$keep * filtered_dt[[exposure]] /
              filtered_dt$prop_scores,
            family = "binomial"
          )
        )

        ## tilt the clever covariate of failure hazard under control
        suppressWarnings(
          fail_haz_0_tilt_fit <- stats::glm(
            filtered_dt$failure ~ -1 + mod_h_0,
            offset = fail_haz_0_est_logit,
            weights = filtered_dt$keep * (1 - filtered_dt[[exposure]]) /
              (1 - filtered_dt$prop_scores),
            family = "binomial"
          )
        )

        ## update the nuisance parameter estimates
        fail_haz_1_est <- predict(fail_haz_1_tilt_fit, type = "response")
        fail_haz_0_est <- predict(fail_haz_0_tilt_fit, type = "response")
        fail_haz_est <- filtered_dt[[exposure]] * fail_haz_1_est +
          (1 - filtered_dt[[exposure]]) * fail_haz_0_est

        ## update the clever covariate
        filtered_dt$failure_haz_est_star <- fail_haz_est
        filtered_dt$failure_haz_exp_est_star <- fail_haz_1_est
        filtered_dt$failure_haz_noexp_est_star <- fail_haz_0_est
        filtered_dt[, `:=`(
          surv_est = cumprod(1 - failure_haz_est_star)
        ), by = "id"]
        filtered_dt[,
                    inner_integrand := keep * ipws / (cens_est_lag * surv_est),
                    by = "id"]
        ## the time cutoff
        filtered_dt[, surv_time_cutoff := min(surv_est)]
        ## compute the partial clever covariates
        filtered_dt[, `:=`(
          partial_h = inner_integrand * surv_time_cutoff,
          partial_h_1 = -surv_time_cutoff * keep /
            (cens_est_lag * surv_est),
          partial_h_0 = surv_time_cutoff * keep /
            (cens_est_lag * surv_est)
        )]
        mod_h <- filtered_dt$partial_h
        mod_h_1 <- filtered_dt$partial_h_1
        mod_h_0 <- filtered_dt$partial_h_0

        ## compute the score
        haz_score <- mod_h * (filtered_dt$failure - fail_haz_est)

        iter <- iter + 1
      }

      # compute the RMST plug in estimator
      filtered_dt[, `:=`(
        surv_exp_est_star = cumprod(1 - failure_haz_exp_est_star),
        surv_noexp_est_star = cumprod(1 - failure_haz_noexp_est_star)
      ), by = "id"]
      filtered_dt <- filtered_dt[filtered_dt[, .I[.N], by = "id"]$V1]
      filtered_dt <- filtered_dt[, `:=`(
        weighted_surv_diff = int_weight *
          (surv_exp_est_star - surv_noexp_est_star)
      )]

      ## return the modifier and surv diff
      filtered_dt[, c("id", "weighted_surv_diff"), with = FALSE]

    }
  )

  ## combine all times and survival probabilities into a long data.table
  survival_preds <- data.table::rbindlist(survival_preds)

  ## compute the difference in restricted mean survival times
  survival_preds[, rmst_diff := sum(weighted_surv_diff), by = "id"]
  estimate <- mean(survival_preds$rmst_diff)

}
