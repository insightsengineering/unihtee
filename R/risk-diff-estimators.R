utils::globalVariables(c(".SD", "..modifiers"))
#' @title One Step Estimator
#'
#' @description \code{one_step_estimator()} computes the one-step estimates of
#'   any given parameter when provided with uncentered efficient influence
#'   functions.
#'
#' @param uncentered_eif_data A \code{data.table} where each column is an
#'   efficient influence function of a potential treatment effect modifier
#'   applied to some dataset.
#'
#' @return A one-row \code{data.table} containing the one-step estimates for
#' each modifier.
#'
#' @importFrom data.table .SD
#'
#' @keywords internal
one_step_estimator <- function(uncentered_eif_data) {
  uncentered_eif_data[, lapply(.SD, mean)]
}

tml_estimator <- function(
  data,
  confounders,
  modifiers,
  exposure,
  outcome,
  prop_score_fit,
  cond_outcome_fit
) {

  # specify the covariates used to predict the potential outcomes
  covariates <- c(confounders, exposure)

  # create dataset where all observations are exposed
  exp_data <- data.table::copy(data)
  exp_data[[exposure]] <- 1

  # create sl3 tasks for exposed dataset
  exp_co_task <- sl3::sl3_Task$new(
    data = exp_data,
    covariates = covariates,
    outcome = outcome
  )
  exp_ps_task <- sl3::sl3_Task$new(
    data = exp_data,
    covariates = confounders,
    outcome = exposure
  )

  # predict outcomes and ps scores of exposed dataset
  exp_outcome <- cond_outcome_fit$fit$predict(exp_co_task)
  exp_ps <- prop_score_fit$fit$predict(exp_ps_task)

  # create dataset where all observations are non-exposed
  noexp_data <- data.table::copy(data)
  noexp_data[[exposure]] <- 0

  # create sl3 tasks for non-exposed dataset
  noexp_co_task <- sl3::sl3_Task$new(
    data = noexp_data,
    covariates = covariates,
    outcome = outcome
  )
  noexp_ps_task <- sl3::sl3_Task$new(
    data = noexp_data,
    covariates = confounders,
    outcome = exposure
  )

  # predict outcomes of non-exposed dataset
  noexp_outcome <- cond_outcome_fit$fit$predict(noexp_co_task)
  noexp_ps <- prop_score_fit$fit$predict(noexp_ps_task)

  # compute that partial clever covariate
  h_partial <- (2 * data[[exposure]] - 1) /
    (data[[exposure]] * prop_score_fit$estimates +
     (1 - data[[exposure]]) * (1 - prop_score_fit$estimates))
  h_partial_1 <- 1 / exp_ps
  h_partial_0 <- -1 / noexp_ps

  # compute TML estimate
  estimates <- lapply(
    modifiers,
    function(mod) {

      # compute the tilted conditional outcome estimators
      mod_var <- var(data[[mod]])
      mod_h <- data[[mod]] * h_partial / mod_var
      mod_h_1 <- data[[mod]] * h_partial_1 / mod_var
      mod_h_0 <- data[[mod]] * h_partial_0 / mod_var
      epsilon <- coef(
        glm(data[[outcome]] ~ -1 + mod_h +
              offset(qlogis(cond_outcome_fit$estimates)),
            family = "quasibinomial")
      )
      q_1_star <- plogis(qlogis(exp_outcome) + epsilon * mod_h_1)
      q_0_star <- plogis(qlogis(noexp_outcome) + epsilon * mod_h_0)

      # compute the plugin estimate with the update cond outcome estimates
      cov(data[[mod]], q_1_star - q_0_star) / mod_var

    })

  # assemble the estimates into a data.table
  names(estimates) <- modifiers
  tmle_dt <- data.table::as.data.table(estimates)
  return(tmle_dt)

}
