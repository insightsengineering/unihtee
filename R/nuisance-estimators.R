#' @title Propensity Score Estimator
#'
#' @description \code{fit_prop_score()} estimates the propensity score over the
#'   \code{train_data} dataset. The estimator used for this estimation is based
#'   on the \code{learners} argument, and the confounders considered are
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
fit_prop_score <- function(
  train_data,
  valid_data,
  learners,
  exposure,
  confounders
) {

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
#'   argument, and the confounders considered are specified by the
#'   \code{confounders} argument.
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
#' @return A named \code{list} of two elements. (1) \code{"estimates"}, the
#'   propensity score estimates for each observation in \code{valid_data}, if
#'   specified, or \code{train_data} otherwise. (2) \code{"fit"}, the trained
#'   \code{\link[sl3]{Stack}} or learner.
#'
#' @importFrom sl3 sl3_Task
#'
#' @keywords internal
fit_cond_outcome <- function(
  train_data,
  valid_data,
  learners,
  outcome,
  exposure,
  confounders
) {

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

    # extract the propensity score predictions for the training data
    cond_outcome_est <- cond_outcome_fit$predict()

  } else {

    # construct the validation data task
    valid_data_task <- sl3::sl3_Task$new(
      data = valid_data,
      covariates = covariates,
      outcome = outcome
    )

   # extract the propensity score predictions for the valid data
    cond_outcome_est <- cond_outcome_fit$predict(valid_data_task)
  }

  return(list(
    "estimates" = cond_outcome_est
  ))

}
