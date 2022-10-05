utils::globalVariables(names = c("a", ".SD"))
#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#' @title
#' @param data
#' @param confounders
#' @param exposure
#' @param outcome
#' @param modifiers
#' @param prop_score_fit
#' @param cond_outcome_fit
#' @return
#' @importFrom data.table copy as.data.table
#' @author boileap2
uncentered_eif <- function(
  data,
  confounders,
  exposure,
  outcome,
  modifiers,
  prop_score_fit,
  prop_score_values,
  cond_outcome_fit
) {

  covariates <- c(confounders, exposure)

  # create dataset where all observations are exposed
  exp_data <- data.table::copy(data)
  exp_data[[exposure]] <- 1

  # create sl3 task for exposed dataset
  exp_task <- sl3::sl3_Task$new(
    data = exp_data,
    covariates = covariates,
    outcome = outcome
  )

  # predict outcomes of exposed dataset
  exp_outcome <- cond_outcome_fit$fit$predict(exp_task)

  # create dataset where all observations are non-exposed
  noexp_data <- data.table::copy(data)
  noexp_data[[exposure]] <- 0

  # create sl3 task for non-exposed dataset
  noexp_task <- sl3::sl3_Task$new(
    data = noexp_data,
    covariates = covariates,
    outcome = outcome
  )

  # predict outcomes of non-exposed dataset
  noexp_outcome <- cond_outcome_fit$fit$predict(noexp_task)

  # compute conditional outcome residuals
  cond_outcome_resid <- data[[outcome]] - cond_outcome_fit$estimates

  # compute the inverse probability weights
  prop_scores <- ifelse(is.null(prop_score_fit),
                        prop_score_values, prop_score_fit$estimates)
  ipws <- (2 * data[[exposure]] - 1) /
    (data[[exposure]] * prop_scores +
      (1 - data[[exposure]]) * (1 - prop_scores))

  # compute augmented inverse probability weights outcomes
  aipws <- ipws * cond_outcome_resid + exp_outcome - noexp_outcome

  # compute that variance of the effect modifiers
  modifier_vars <- sapply(modifiers, function(modifier) var(data[[modifier]]))

  # vectorized uncentered eifs
  modifiers_dt <- lapply(
    modifiers,
    function(modifier) {
      mod_vec <- data[[modifier]]
      mod_vec - mean(mod_vec)
  })
  modifiers_dt <- data.table::as.data.table(modifiers_dt)
  colnames(modifiers_dt) <- modifiers
  eif_dt <- tcrossprod(aipws, 1 / modifier_vars) * modifiers_dt

  return(eif_dt)

}
