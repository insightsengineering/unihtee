utils::globalVariables(
  c(".SD", "..modifiers", ".I", ".N", "ipws", "surv_est", "cens_est",
    "cens_est_lag", "surv_exp_est", "surv_noexp_est", "surv_exp_est_star",
    "surv_noexp_est_star", "prev_time", "int_weight", "inner_integrand",
    "surv_time_cutoff", "partial_h", "partial_h_1", "partial_h_0",
    "failure_haz_est_star", "failure_haz_exp_est", "failure_haz_noexp_est",
    "failure_haz_exp_est_star", "failure_haz_noexp_est_star")
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
#' @param effect A \code{character} indicating the type of treatment effect
#'   modifier variable importance parameter. Currently supports
#'   \code{"absolute"} and \code{"relative"}.
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
                          effect,
                          prop_score_fit,
                          prop_score_values = NULL,
                          cond_outcome_fit,
                          failure_hazard_fit,
                          censoring_hazard_fit) {

  # constant used to truncate estimates near zero
  eps <- .Machine$double.eps

  ## compute the inverse probability weights
  ## NOTE: PS estimates must be as long as long_dt in TTE settings
  if (!is.null(prop_score_values)) {
    prop_scores <- data[[prop_score_values]]
  } else {
    prop_scores <- prop_score_fit$estimates
  }

  if (!is.null(cond_outcome_fit)) {
    # truncate cond outcome estimates if necessary
    estimates <- cond_outcome_fit$estimates
    exp_estimates <- cond_outcome_fit$exp_estimates
    noexp_estimates <- cond_outcome_fit$noexp_estimates
    if (effect == "relative") {
      estimates[estimates < eps] <- eps
      exp_estimates[exp_estimates < eps] <- eps
      noexp_estimates[noexp_estimates < eps] <- eps
    }

    # compute that partial clever covariate
    h_partial_1 <- data[[exposure]] / prop_scores
    h_partial_0 <- (1 - data[[exposure]]) / (1 - prop_scores)

    ## set the TML tilting hyperparameters
    score_stop_crit <- 1 / nrow(data) * log(nrow(data))
    tilt_tol <- 5
    max_iter <- 10

    ## compute TML estimate
    estimates <- lapply(
      modifiers,
      function(mod) {

        ## compute the tilted conditional outcome estimators
        centered_mod <- scale(data[[mod]], scale = FALSE, center = TRUE)
        attr(centered_mod, "scaled:center") <- NULL
        centered_mod <- as.vector(centered_mod)
        mod_var <- var(centered_mod)

        ## tilt the estimators
        q_score <- Inf
        iter <- 1
        q_star <- estimates
        q_1_star <- exp_estimates
        q_0_star <- noexp_estimates
        if (effect == "absolute") {
          scaled_mod <- centered_mod / mod_var
        } else if (effect == "relative") {
          scaled_mod <- centered_mod / (mod_var * q_star)
        }

        while (score_stop_crit < abs(mean(q_score)) && iter < max_iter) {

          ## transform expected conditional outcomes for tilting
          q_star_logit <- stats::qlogis(bound_precision(q_star))
          q_1_star_logit <- stats::qlogis(bound_precision(q_1_star))
          q_0_star_logit <- stats::qlogis(bound_precision(q_0_star))

          ## tilt the clever covariate using the log loss for q_1
          suppressWarnings(
            q_1_tilt_fit <- stats::glm(
              data[[outcome]] ~ -1 + scaled_mod,
              offset = q_star_logit,
              weights = h_partial_1,
              family = "binomial",
              start = 0
            )
          )
          if (is.na(stats::coef(q_1_tilt_fit))) {
            q_1_tilt_fit$coefficients <- 0
          } else if (!q_1_tilt_fit$converged ||
                       abs(max(stats::coef(q_1_tilt_fit))) > tilt_tol) {
            q_1_tilt_fit$coefficients <- 0
          }
          epsilon_1 <- unname(stats::coef(q_1_tilt_fit))

          ## tilt the clever covariate using the log loss for q_0
          suppressWarnings(
            q_0_tilt_fit <- stats::glm(
              data[[outcome]] ~ -1 + scaled_mod,
              offset = q_star_logit,
              weights = h_partial_0,
              family = "binomial",
              start = 0
            )
          )
          if (is.na(stats::coef(q_0_tilt_fit))) {
            q_0_tilt_fit$coefficients <- 0
          } else if (!q_0_tilt_fit$converged ||
                       abs(max(stats::coef(q_0_tilt_fit))) > tilt_tol) {
            q_0_tilt_fit$coefficients <- 0
          }
          epsilon_0 <- unname(stats::coef(q_0_tilt_fit))

          ## update the nuisance parameter estimates
          q_1_star <- stats::plogis(q_1_star_logit + epsilon_1 * scaled_mod)
          q_0_star <- stats::plogis(q_0_star_logit + epsilon_0 * scaled_mod)
          q_star <- data[[exposure]] * q_1_star +
            (1 - data[[exposure]]) * q_0_star

          ## update the relative TEM-VIP's clever covariate
          if (effect == "relative") {
            ## bound q, just in case
            q_star <- bound_precision(q_star)
            q_1_star <- bound_precision(q_1_star)
            q_0_star <- bound_precision(q_0_star)
            scaled_mod <- centered_mod / (mod_var * q_star)
          }

          ## compute the score
          if (effect == "absolute") {
            q_score <- centered_mod * (2 * data[[exposure]] - 1) /
              (data[[exposure]] * prop_scores +
                 (1 - data[[exposure]]) * (1 - prop_scores)) /
              mod_var * (data[[outcome]] - q_star)
          } else if (effect == "relative") {
            q_score <- centered_mod * (2 * data[[exposure]] - 1) /
              (data[[exposure]] * prop_scores +
                 (1 - data[[exposure]]) * (1 - prop_scores)) /
              mod_var * (data[[outcome]] - q_star) / q_star
          }

          iter <- iter + 1

        }

        ## compute the plugin estimate with the update cond outcome estimates
        if (effect == "absolute") {
          stats::cov(centered_mod, q_1_star - q_0_star) / mod_var
        } else if (effect == "relative") {
          stats::cov(centered_mod, log(q_1_star) - log(q_0_star)) / mod_var
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

    ## absolute TEM VIP
    if (effect == "absolute") {
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

          tilted_survivals <- lapply(
            modifiers,
            function(mod) {

              ## center the modifier
              centered_mod <- scale(
                filtered_dt[[mod]], scale = FALSE, center = TRUE
              )
              attr(centered_mod, "scaled:center") <- NULL
              centered_mod <- as.vector(centered_mod)

              ## finalize the clever covariate
              mod_h <- centered_mod * filtered_dt$partial_h /
                mod_vars[[mod]]
              mod_h_1 <- centered_mod * filtered_dt$partial_h_1 /
                mod_vars[[mod]]
              mod_h_0 <- centered_mod * filtered_dt$partial_h_0 /
                mod_vars[[mod]]

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
                    offset = fail_haz_est_logit,
                    weights = filtered_dt$keep * filtered_dt[[exposure]] /
                      filtered_dt$prop_scores,
                    family = "binomial"
                  )
                 )
                if (is.na(stats::coef(fail_haz_1_tilt_fit))) {
                  fail_haz_1_tilt_fit$coefficients <- 0
                } else if (!fail_haz_1_tilt_fit$converged ||
                            abs(max(stats::coef(fail_haz_1_tilt_fit)))
                           > tilt_tol) {
                  fail_haz_1_tilt_fit$coefficients <- 0
                }
                epsilon_1 <- unname(stats::coef(fail_haz_1_tilt_fit))

                ## tilt the clever covariate of failure hazard under control
                suppressWarnings(
                 fail_haz_0_tilt_fit <- stats::glm(
                    filtered_dt$failure ~ -1 + mod_h_0,
                    offset = fail_haz_est_logit,
                    weights = filtered_dt$keep * (1 - filtered_dt[[exposure]]) /
                      (1 - filtered_dt$prop_scores),
                    family = "binomial"
                  )
                 )
                if (is.na(stats::coef(fail_haz_0_tilt_fit))) {
                  fail_haz_0_tilt_fit$coefficients <- 0
                } else if (!fail_haz_0_tilt_fit$converged ||
                            abs(max(stats::coef(fail_haz_0_tilt_fit)))
                           > tilt_tol) {
                  fail_haz_0_tilt_fit$coefficients <- 0
                }
                epsilon_0 <- unname(stats::coef(fail_haz_0_tilt_fit))


                ## update the nuisance parameter estimates
                fail_haz_1_est <- bound_precision(
                  stats::plogis(fail_haz_1_est_logit + epsilon_1 * mod_h_1)
                )
                fail_haz_0_est <- bound_precision(
                  stats::plogis(fail_haz_0_est_logit + epsilon_0 * mod_h_0)
                )
                fail_haz_est <- filtered_dt[[exposure]] * fail_haz_1_est +
                  (1 - filtered_dt[[exposure]]) * fail_haz_0_est

                ## update the clever covariate
                filtered_dt$failure_haz_est_star <- fail_haz_est
                filtered_dt$failure_haz_exp_est_star <- fail_haz_1_est
                filtered_dt$failure_haz_noexp_est_star <- fail_haz_0_est
                filtered_dt[, `:=`(
                  surv_est = cumprod(1 - failure_haz_est_star),
                  surv_exp_est_star = cumprod(1 - failure_haz_exp_est_star),
                  surv_noexp_est_star = cumprod(1 - failure_haz_noexp_est_star)
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
                mod_h <- centered_mod * filtered_dt$partial_h /
                  mod_vars[[mod]]
                mod_h_1 <- centered_mod * filtered_dt$partial_h_1 /
                  mod_vars[[mod]]
                mod_h_0 <- centered_mod * filtered_dt$partial_h_0 /
                  mod_vars[[mod]]

                ## compute the score
                haz_score <- mod_h * (filtered_dt$failure - fail_haz_est)

                iter <- iter + 1
              }

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

    } else if (effect == "relative") {
      ## compute the clever covariate
      data[, `:=`(
        partial_h = keep * ipws / (cens_est_lag * surv_est),
        partial_h_1 = keep / (cens_est_lag * surv_est * prop_scores),
        partial_h_0 = -keep / (cens_est_lag * surv_est * (1 - prop_scores))
      )]

      estimates <- lapply(
        modifiers,
        function(mod) {

          data_mod <- data.table::copy(data)

          ## center the modifier
          centered_mod <- scale(data_mod[[mod]], scale = FALSE, center = TRUE)
          attr(centered_mod, "scaled:center") <- NULL
          centered_mod <- as.vector(centered_mod)

          ## finalize the clever covariate
          mod_h <- centered_mod * data_mod$partial_h / mod_vars[[mod]]
          mod_h_1 <- centered_mod * data_mod$partial_h_1 / mod_vars[[mod]]
          mod_h_0 <- centered_mod * data_mod$partial_h_0 / mod_vars[[mod]]

          ## set initial values
          data_mod$failure_haz_est_star <- data_mod$failure_haz_est
          data_mod$failure_haz_exp_est_star <- data_mod$failure_haz_exp_est
          data_mod$failure_haz_noexp_est_star <- data_mod$failure_haz_noexp_est
          epsilon <- 1
          iter <- 1
          max_iter <- 5
          haz_score <- Inf
          score_stop_crit <- 1 / (nrow(data_mod) * log(nrow(data_mod)))

          ## bound q, just in case
          data_mod$failure_haz_est_star[
            data_mod$failure_haz_est_star < (0 + eps)
          ] <- 0 + eps
          data_mod$failure_haz_est_star[
            data_mod$failure_haz_est_star > (1 - eps)
          ] <- 1 - eps
          data_mod$failure_haz_exp_est_star[
            data_mod$failure_haz_exp_est_star < (0 + eps)
          ] <- 0 + eps
          data_mod$failure_haz_exp_est_star[
            data_mod$failure_haz_exp_est_star > (1 - eps)
          ] <- 1 - eps
          data_mod$failure_haz_noexp_est_star[
            data_mod$failure_haz_noexp_est_star < (0 + eps)
          ] <- 0 + eps
          data_mod$failure_haz_noexp_est_star[
            data_mod$failure_haz_noexp_est_star > (1 - eps)
          ] <- 1 - eps

          ## tilt the conditional failure estimates
          while (abs(mean(haz_score)) > score_stop_crit && iter < max_iter) {
            epsilon <- stats::coef(
              stats::glm(
                data_mod$failure ~ -1 + mod_h,
                offset = data_mod$failure_haz_est_star,
                family = "gaussian"
              )
            )
            ## update hazard
            data_mod$failure_haz_est_star <- data_mod$failure_haz_est +
              epsilon * mod_h

            ## compute the tilted conditional survival probabilities differences
            data_mod$failure_haz_exp_est_star <-
              data_mod$failure_haz_exp_est_star + epsilon * mod_h_1
            data_mod$failure_haz_noexp_est_star <-
              data_mod$failure_haz_noexp_est_star + epsilon * mod_h_0

            ## compute the score
            haz_score <- mod_h *
              (data_mod$failure - data_mod$failure_haz_est_star)

            ## bound q, just in case
            data_mod$failure_haz_est_star[
              data_mod$failure_haz_est_star < (0 + eps)] <- 0 + eps
            data_mod$failure_haz_est_star[
              data_mod$failure_haz_est_star > (1 - eps)] <- 1 - eps
            data_mod$failure_haz_exp_est_star[
              data_mod$failure_haz_exp_est_star < (0 + eps)] <- 0 + eps
            data_mod$failure_haz_exp_est_star[
              data_mod$failure_haz_exp_est_star > (1 - eps)] <- 1 - eps
            data_mod$failure_haz_noexp_est_star[
              data_mod$failure_haz_noexp_est_star < (0 + eps)] <- 0 + eps
            data_mod$failure_haz_noexp_est_star[
              data_mod$failure_haz_noexp_est_star > (1 - eps)] <- 1 - eps

            iter <- iter + 1
          }

          ## compute the conditional survival under treatment and control
          data_mod[, `:=`(
            surv_exp_est_star = cumprod(1 - failure_haz_exp_est_star),
            surv_noexp_est_star = cumprod(1 - failure_haz_noexp_est_star)
          ), by = "id"]

          ## only retain the rows at the designated time cutoff
          data_mod <- data_mod[data_mod[, .I[.N], by = "id"]$V1]

          ## compute the estimate
          log_diff <- log(data_mod$surv_exp_est_star /
                            data_mod$surv_noexp_est_star)
          cov(log_diff, data_mod[[mod]]) / mod_vars[[mod]]

        }
      )
    }
  }

  # assemble the estimates into a data.table
  names(estimates) <- modifiers
  tmle_dt <- data.table::as.data.table(estimates)
  return(tmle_dt)
}


###############################################################################

#' Bounding to numerical precision
#'
#' Bounds extreme values to a specified tolerance level, for use with sensitive
#' quantities that must be transformed, e.g., via \code{\link[stats]{qlogis}}.
#'
#' @param vals A \code{numeric} vector of values in the unit interval \[0, 1\].
#' @param tol A \code{numeric} indicating the tolerance limit to which extreme
#'  values should be truncated. Realizations of \code{val} less than \code{tol}
#'  are truncated to \code{tol} while those greater than (1 - \code{tol}) are
#'  truncated to (1 - \code{tol}).
#'
#' @importFrom assertthat assert_that
#'
#' @author Nima S. Hejazi
#'
#' @keywords internal
bound_precision <- function(vals, tol = 1e-6) {
  vals[vals < tol] <- tol
  vals[vals > 1 - tol] <- 1 - tol
  return(vals)
}
