utils::globalVariables(names = c("a", ".SD", "surv_est"))
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
#' @param type A \code{character} indicating the type of treatment effect
#'   modifier variable importance parameter. Currently supports
#'   \code{"risk difference"} and \code{"relative risk"}.
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
#'   \code{\link{fit_failure_hazard}()} function.'
#' @param failure_hazard_fit A \code{list} output by the
#'   \code{\link{fit_cond_outcome}()} function.'
#' @param censoring_hazard_fit A \code{list} output by the
#'   \code{\link{fit_censoring_hazard}()} function.'
#'
#' @return A \code{data.table} whose columns are the uncentered efficient
#'   influence functions of each variable in \code{modifiers}. The rows
#'   correspond to the observations of \code{data}.
#'
#' @importFrom data.table copy as.data.table
#'
#' @keywords internal

uncentered_eif <- function(data,
                           type,
                           confounders,
                           exposure,
                           outcome,
                           modifiers,
                           prop_score_fit,
                           prop_score_values,
                           cond_outcome_fit,
                           failure_hazard_fit,
                           censoring_hazard_fit) {

  ## compute the inverse probability weights
  ## NOTE: PS estimates must be as long as long_dt in TTE settings
  if (!is.null(prop_score_values)) {
    prop_scores <- prop_score_values
  } else {
    prop_scores <- prop_score_fit$estimates
  }
  ipws <- (2 * data[[exposure]] - 1) /
    (data[[exposure]] * prop_scores +
      (1 - data[[exposure]]) * (1 - prop_scores))

  ## decide which EIFs to compute
  ## EIFs for continuous and binary outcomes with binary exposure
  if (is.null(failure_hazard_fit) && is.null(censoring_hazard_fit)) {

    ## compute conditional outcome residuals
    cond_outcome_resid <- data[[outcome]] - cond_outcome_fit$estimates

    ## compute augmented inverse probability weights outcomes
    if (type == "risk difference") {
      aipws <- ipws * cond_outcome_resid + cond_outcome_fit$exp_estimates -
        cond_outcome_fit$noexp_estimates
    } else if (type == "relative risk") {
      ## NOTE: Make sure not to divide by zero or compute log of zero
      eps <- 1e-10
      estimates <- cond_outcome_fit$estimates
      estimates[estimates < eps] <- eps
      exp_estimates <- cond_outcome_fit$exp_estimates
      exp_estimates[exp_estimates < eps] <- eps
      noexp_estimates <- cond_outcome_fit$noexp_estimates
      noexp_estimates[noexp_estimates < eps] <- eps
      aipws <- ipws * cond_outcome_resid / estimates +
        log(exp_estimates) - log(noexp_estimates)
    }

    ## EIFs for time-to-event outcomes
  } else if (is.null(cond_outcome_fit)) {
    ## compute the survival probabilities
    data$ipws <- ipws
    data$failure_haz_est <- failure_hazard_fit$estimates
    data$failure_haz_exp_est <- failure_hazard_fit$exp_estimates
    data$failure_haz_noexp_est <- failure_hazard_fit$noexp_estimates
    data$censoring_haz_est <- censoring_hazard_fit$estimates
    data[, surv_est := cumprod(1 - failure_haz_est), by = "id"]
    data[, `:=`(
      surv_time_cutoff = min(surv_est),
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
    data[, int_exp := int_weight * keep * ipws * surv_time_cutoff /
           (cens_est_lag * surv_est) * (failure - failure_haz_est),
         by = "id"]
    data[, int_exp := cumsum(int_exp), by = "id"]

    ## use only the observations at the time_cutoff
    time_cutoff <- max(data[[outcome]])
    data <- data[get(outcome) == time_cutoff, ]

    ## compute the uncentered EIF
    aipws <- data$int_exp + data$surv_exp_est - data$surv_noexp_est

  }

  ## compute that variance of the effect modifiers
  modifier_vars <- sapply(modifiers, function(modifier) var(data[[modifier]]))

  ## vectorized uncentered eifs
  modifiers_dt <- lapply(
    modifiers,
    function(modifier) {
      mod_vec <- data[[modifier]]
      mod_vec - mean(mod_vec)
    }
  )
  modifiers_dt <- data.table::as.data.table(modifiers_dt)
  colnames(modifiers_dt) <- modifiers
  eif_dt <- tcrossprod(aipws, 1 / modifier_vars) * modifiers_dt

  return(eif_dt)
}
