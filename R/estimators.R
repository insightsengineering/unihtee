utils::globalVariables(
  c(".SD", "..modifiers", ".I", ".N", "ipws", "surv_est", "cens_est",
    "cens_est_lag", "surv_exp_est", "surv_noexp_est", "surv_exp_est_star",
    "surv_noexp_est_star", "prev_time", "int_weight", "inner_integrand",
    "surv_time_cutoff", "partial_h", "partial_h_1", "partial_h_0",
    "failure_haz_exp_est", "failure_haz_noexp_est", "failure_haz_exp_est_star",
    "failure_haz_noexp_est_star")
)
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
#' each potential modifier.
#'
#' @importFrom data.table .SD
#'
#' @keywords internal
one_step_estimator <- function(uncentered_eif_data) {
  uncentered_eif_data[, lapply(.SD, mean)]
}

#' @title Targeted Maximum Likelihood Estimator
#'
#' @description \code{tml_estimator()} computes the targeted maximum likelihood
#'   estimates of any given parameter.
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
#' @param type A \code{character} indicating the type of treatment effect
#'   modifier variable importance parameter. Currently supports
#'   \code{"risk difference"} and \code{"relative risk"}.
#' @param prop_score_fit A \code{list} output by the
#'   \code{\link{fit_prop_score}()} function.
#' @param prop_score_values A \code{numeric} vector corresponding to the (known)
#'   propensity score values for each observation in \code{data}.
#' @param cond_outcome_fit A \code{list} output by the
#'   \code{\link{fit_cond_outcome}()} function.
#' @param failure_hazard_fit A \code{list} output by the
#'   \code{\link{fit_cond_outcome}()} function.'
#' @param censoring_hazard_fit A \code{list} output by the
#'   \code{\link{fit_censoring_hazard}()} function.'
#'
#' @return A one-row \code{data.table} containing the targeted maximum
#'   likelihood estimates for each potential modifier.
#'
#' @importFrom data.table as.data.table .I
#' @importFrom stats glm qlogis plogis coef cov
#'
#' @keywords internal
tml_estimator <- function(data,
                          confounders,
                          modifiers,
                          exposure,
                          outcome,
                          type,
                          prop_score_fit,
                          prop_score_values = NULL,
                          cond_outcome_fit,
                          failure_hazard_fit,
                          censoring_hazard_fit) {

  # constant used to truncate estimates near zero
  eps <- 1e-10

  ## compute the inverse probability weights
  ## NOTE: PS estimates must be as long as long_dt in TTE settings
  if (!is.null(prop_score_values)) {
    prop_scores <- prop_score_values
  } else {
    prop_scores <- prop_score_fit$estimates
  }

  if (!is.null(cond_outcome_fit)) {
    # truncate cond outcome estimates if necessary
    estimates <- cond_outcome_fit$estimates
    exp_estimates <- cond_outcome_fit$exp_estimates
    noexp_estimates <- cond_outcome_fit$noexp_estimates
    if (type == "relative risk") {
      estimates[estimates < eps] <- eps
      exp_estimates[exp_estimates < eps] <- eps
      noexp_estimates[noexp_estimates < eps] <- eps
    }

    # compute that partial clever covariate
    if (type == "risk difference") {
      h_partial <- (2 * data[[exposure]] - 1) /
        (data[[exposure]] * prop_scores +
         (1 - data[[exposure]]) * (1 - prop_scores))
      h_partial_1 <- 1 / prop_scores
      h_partial_0 <- -1 / (1 - prop_scores)
    } else if (type == "relative risk") {
      # NOTE: Make sure not to divide by zero... or take log of zero
      h_partial <- (2 * data[[exposure]] - 1) /
        ((data[[exposure]] * prop_scores +
          (1 - data[[exposure]]) * (1 - prop_scores)) * estimates)
      h_partial_1 <- 1 / (prop_scores * exp_estimates)
      h_partial_0 <- -1 / ((1 - prop_scores) * noexp_estimates)
    }

    ## compute TML estimate
    estimates <- lapply(
      modifiers,
      function(mod) {

        ## compute the tilted conditional outcome estimators
        mod_var <- var(data[[mod]])
        mod_h <- data[[mod]] * h_partial / mod_var
        mod_h_1 <- data[[mod]] * h_partial_1 / mod_var
        mod_h_0 <- data[[mod]] * h_partial_0 / mod_var

        ## tilt the estimators
        epsilon <- 1
        q_star <- estimates
        q_1_star <- exp_estimates
        q_0_star <- noexp_estimates
        while(abs(epsilon) > 1e-10) {
          epsilon <- stats::coef(
            stats::glm(
              data[[outcome]] ~ -1 + mod_h,
              offset = stats::qlogis(q_star),
              family = "quasibinomial"
            )
          )
          q_star <- stats::plogis(
            stats::qlogis(q_star) + epsilon * mod_h
          )
          q_1_star <- stats::plogis(
            stats::qlogis(q_1_star) + epsilon * mod_h_1
          )
          q_0_star <- stats::plogis(
            stats::qlogis(q_0_star) + epsilon * mod_h_0
          )
        }

        ## compute the plugin estimate with the update cond outcome estimates
        if (type == "risk difference") {
          stats::cov(data[[mod]], q_1_star - q_0_star) / mod_var
        } else if (type == "relative risk") {
          stats::cov(data[[mod]], log(q_1_star) - log(q_0_star)) / mod_var
        }
      }
    )

  } else if (is.null(cond_outcome_fit)) {

    ## compute modifier variances
    mod_dt <- data[data[, .I[1], by = "id"]$V1]
    mod_dt <- mod_dt[, ..modifiers]
    mod_vars <- mod_dt[, lapply(.SD, var)]

    ## compute the IPWs
    data$prop_scores <- prop_scores
    data$ipws <- (2 * data[[exposure]] - 1) /
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

    # compute the partial clever covariates at each timepoint
    the_times <- unique(data[[outcome]])
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
          partial_h_1 = surv_time_cutoff * keep /
            (cens_est_lag * surv_est * prop_scores),
          partial_h_0 = -surv_time_cutoff * keep /
            (cens_est_lag * surv_est * (1 - prop_scores))
        )]

        tilted_survivals <- lapply(
          modifiers,
          function(mod) {

            ## finalize the clever covariate
            mod_h <- filtered_dt[[mod]] * filtered_dt$partial_h /
              mod_vars[[mod]]
            mod_h_1 <- filtered_dt[[mod]] * filtered_dt$partial_h_1 /
              mod_vars[[mod]]
            mod_h_0 <- filtered_dt[[mod]] * filtered_dt$partial_h_0 /
              mod_vars[[mod]]

            ## set initial values
            init_fail_haz_est <- filtered_dt$failure_haz_est
            filtered_dt$failure_haz_exp_est_star <-
              filtered_dt$failure_haz_exp_est
            filtered_dt$failure_haz_noexp_est_star <-
              filtered_dt$failure_haz_noexp_est
            epsilon <- 1

            ## tilt the conditional failure estimates
            while (abs(epsilon) > 1e-10) {
              epsilon <- stats::coef(
                stats::glm(
                  filtered_dt$failure ~ -1 + mod_h,
                  offset = stats::qlogis(init_fail_haz_est),
                  family = "quasibinomial"
                )
              )
              init_fail_haz_est <- stats::plogis(
                stats::qlogis(init_fail_haz_est) + epsilon * mod_h
              )

              ## compute the tilted conditional survival probabilities
              ## differences
              filtered_dt$failure_haz_exp_est_star <- stats::plogis(
                stats::qlogis(filtered_dt$failure_haz_exp_est_star) +
                  epsilon * mod_h_1
              )
              filtered_dt$failure_haz_noexp_est_star <- stats::plogis(
                stats::qlogis(filtered_dt$failure_haz_noexp_est_star) +
                  epsilon * mod_h_0
              )
            }

            filtered_dt[, `:=`(
              surv_exp_est_star = cumprod(1 - failure_haz_exp_est_star),
              surv_noexp_est_star = cumprod(1 - failure_haz_noexp_est_star)
            ), by = "id"]
            filtered_dt <- filtered_dt[filtered_dt[, .I[.N], by = "id"]$V1]
            filtered_dt <- filtered_dt[, `:=`(
              weighted_surv_diff = int_weight *
                (surv_exp_est_star - surv_noexp_est_star),
              modifier = mod
            )]

            ## return the modifier and surv diff
            filtered_dt[, c("id", "modifier", "weighted_surv_diff"),
                        with = FALSE]

          }
        )

        ## combine all the modifiers into a long data.table
        data.table::rbindlist(tilted_survivals)

      }
    )

    ## combine all times and survival probabilities into a long data.table
    survival_preds <- data.table::rbindlist(survival_preds)

    ## compute the restricted mean survival times
    rmst_dt <- data.table::dcast(
      survival_preds, id ~ modifier, value.var = "weighted_surv_diff",
      fun.aggregate = sum
    )

    ## compute the TML estimates
    estimates <- lapply(
      modifiers,
      function(mod) {
        cov(rmst_dt[[mod]], mod_dt[[mod]]) / mod_vars[[mod]]
      }
    )

  }

  # assemble the estimates into a data.table
  names(estimates) <- modifiers
  tmle_dt <- data.table::as.data.table(estimates)
  return(tmle_dt)
}
