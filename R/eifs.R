utils::globalVariables(names = c("a", ".SD"))
#' @title Uncentered Efficient Influence Function Computer
#'
#' @description \code{uncentered_eif()} computes thei efficient influence
#'   function for the chosen parameter using the already estimated nuisance
#'   parameters. If certain nuisance parameters are known, such as propensity
#'   scores in a randomized control trial, then they may be input directly into
#'   this function.
#'
#' @param data A \code{data.table} containing the observed data.
#'   \code{train_data} is formatted by \code{\link{unihtee}()}.
#' @param confounders A \code{character} vector of column names corresponding to
#'   baseline covariates.
#' @param exposure A \code{character} corresponding to the exposure variable.
#' @param outcome A \code{character} corresponding to the outcome variable.
#' @param modifiers A \code{character} vector of columns names corresponding to
#'   the suspected effect modifiers. This vector must be a subset of
#'   \code{confounders}.
#' @param prop_score_fit A \code{list} output by the
#'   \code{\link{fit_prop_score}()} function.
#' @param prop_score_values A \code{numeric} vector corresponding to the (known)
#'   propensity score values for each observation in \code{data}.
#' @param cond_outcome_fit A \code{list} output by the
#'   \code{\link{fit_cond_outcome}()} function.
#'
#' @return A \code{data.table} whose columns are the uncentered efficient
#'   influence functions of each variable in \code{modifiers}. The rows
#'   correspond to the observations of \code{data}.
#'
#' @importFrom data.table copy as.data.table
#'
#' @keywords internal

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

  # specify the covariates used to predict the potential outcomes
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
  prop_scores <- ifelse(!is.null(prop_score_values),
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
