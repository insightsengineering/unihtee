utils::globalVariables(c("..to_keep", ".SD", ".I", "p_value"))
#' @title Univariate Heterogeneous Treatment Effect Modifier Estimator
#'
#' @description \code{unihtee()} estimates treatment effect modifiers variable
#'   importance parameters (TEM VIPs). These TEM VIPs are defined on the
#'   absolute or relative scales and can be estimated using one-step or targeted
#'   maximum likelihood estimators.
#'
#' @param data A \code{data.table} containing the observed data.
#' @param confounders A \code{character} vector of column names corresponding to
#'   baseline covariates.
#' @param modifiers A \code{character} vector of columns names corresponding to
#'   the suspected effect modifiers. This vector must be a subset of
#'   \code{confounders}.
#' @param exposure A \code{character} corresponding to the exposure variable.
#' @param outcome A \code{character} corresponding to the outcome variable.
#' @param censoring A \code{character} indicating the right censoring indicator
#'   variable. Only used with time-to-event outcomes. Defaults to \code{NULL}.
#' @param time_cutoff A \code{numeric} representing the timepoint at which to
#'   evaluate the time-to-event parameter. Only used with time-to-event
#'   outcomes. Defaults to \code{NULL}.
#' @param outcome_type A \code{character} indicating the outcome type.
#'   \code{"continuous"}, \code{"binary"} and \code{"time-to-event"} are
#'   currently supported.
#' @param scale A \code{character} indicating the type of treatment effect
#'   modifier variable importance parameter. Currently supports
#'   \code{"absolute"} and \code{"relative"}.
#' @param estimator A \code{character} set to either \code{"tmle"} or
#'   \code{"onestep"}. The former results in \code{unihtee()} to use a targeted
#'   maximum likelihood estimators to estimate the deisred TEM VIP, while the
#'   latter uses a one step estimator.
#' @param cross_fit A \code{logical} determining whether cross-fitting should be
#'   used. Defaults to \code{FALSE}.
#' @param cross_fit_folds A \code{numeric} stating the number of folds to use in
#'   the cross-fitting procedure. Defaults to 5.
#' @param cond_outcome_estimator A \code{\link[sl3]{Stack}}, or other learner
#'   class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
#'   learners from \pkg{sl3} to estimate the conditional outcome. Defaults to a
#'   generalized linear model with one- and two- way interactions among all
#'   \code{confounders} and \code{exposure} variables. Only used with continuous
#'   and binary outcomes.
#' @param prop_score_estimator A \code{\link[sl3]{Stack}}, or other learner
#'   class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
#'   learners from \pkg{sl3} to estimate the propensity score. Defaults to a
#'   generalized linear model with one- and two- way interactions among all
#'   \code{confounders} variables.
#' @param prop_score_values An optional \code{character} corresponding to the
#'   (known) propensity score values for each observation in \code{data}.
#'   Defaults to \code{NULL}.
#' @param failure_hazard_estimator A \code{\link[sl3]{Stack}}, or other learner
#'   class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
#'   learners from \pkg{sl3} to estimate the conditional failure hazard
#'   function. Defaults to an xgboost learner with \code{confounders} and
#'   \code{exposure} variables as covariates. Only used with time-to-event
#'   outcomes.
#' @param censoring_hazard_estimator A \code{\link[sl3]{Stack}}, or other
#'   learner class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a
#'   set of learners from \pkg{sl3} to estimate the conditional censoring hazard
#'   function. Defaults to an xgboost learner with \code{confounders} and
#'   \code{exposure} variables as covariates. Only used with time-to-event
#'   outcomes.
#' @param parallel A \code{logical} stating if
#'   \code{\link[origami:cross_validate]{origami}}'s built-in parallelized
#'   cross-validation routines should be used when \code{cross_fit = TRUE}. The
#'   \href{https://cran.r-project.org/package=future}{\code{future}} suite is
#'   used. Defaults to \code{FALSE}.
#'
#'
#' @return A \code{data.table} containing the effect estimates and (adjusted)
#'   p-values of the \code{modifiers}. The suspected treatment effect modifiers
#'   ordered according to ascending p-values.
#'
#' @importFrom data.table as.data.table .I .SD rbindlist
#'
#' @export
unihtee <- function(data,
                    confounders,
                    modifiers,
                    exposure,
                    outcome,
                    censoring = NULL,
                    time_cutoff = NULL,
                    outcome_type = c("continuous", "binary", "time-to-event"),
                    scale = c("absolute", "relative"),
                    estimator = c("tmle", "onestep"),
                    cross_fit = FALSE,
                    cross_fit_folds = 5,
                    cond_outcome_estimator = sl3::Lrnr_glm_fast$new(),
                    prop_score_estimator = sl3::Lrnr_glm_fast$new(),
                    prop_score_values = NULL,
                    failure_hazard_estimator = sl3::Lrnr_xgboost$new(),
                    censoring_hazard_estimator = sl3::Lrnr_xgboost$new(),
                    parallel = FALSE) {

  ## specify the TEM VIP type
  param_scale <- match.arg(scale)

  ## specify the outcome type
  outcome_type <- match.arg(outcome_type)

  ## set the estimator
  estimator <- match.arg(estimator)

  ## transform data into data.table object with required columns
  data <- data.table::as.data.table(data)
  to_keep <- unique(c(
    confounders, modifiers, exposure, outcome, prop_score_values, censoring
  ))
  data <- data[, ..to_keep]

  ## get the number of observations
  n_obs <- nrow(data)

  ## scale the outcome to be between 0 and 1 if outcome is continuous
  ## NOTE: Can remove this if using linear update
  if (outcome_type == "never") {
    min_out <- min(data[[outcome]])
    max_out <- max(data[[outcome]])
    rescale_factor <- max_out - min_out
    data[[outcome]] <- (data[[outcome]] - min_out) / rescale_factor
  } else {
    rescale_factor <- 1
  }

  ## create the long-form data for time-to-event outcomes
  if (outcome_type == "time-to-event") {
    data <- tte_data_melt(
      data = data,
      confounders = confounders,
      exposure = exposure,
      outcome = outcome,
      censoring = censoring,
      time_cutoff = time_cutoff,
      prop_score_values = prop_score_values
    )
  }

  if (!cross_fit) {

    ## fit the propensity score, if necessary
    if (outcome_type == "time-to-event") {
      prop_score_train_data <- data[data[, .I[1], by = "id"]$V1]
    } else {
      prop_score_train_data <- data.table::copy(data)
    }
    if (is.null(prop_score_values)) {
      prop_score_fit <- fit_prop_score(
        train_data = prop_score_train_data,
        valid_data = data,
        learners = prop_score_estimator,
        exposure = exposure,
        confounders = confounders
      )
    } else {
      prop_score_fit <- list(
        "estimates" = data[[prop_score_values]]
      )
    }

    ## fit the remaining nuisance parameter for continuous and binary outcomes
    if (outcome_type %in% c("continuous", "binary")) {
      cond_outcome_fit <- fit_cond_outcome(
        train_data = data,
        valid_data = NULL,
        learners = cond_outcome_estimator,
        exposure = exposure,
        confounders = confounders,
        outcome = outcome
      )
      failure_hazard_fit <- NULL
      censoring_hazard_fit <- NULL
    }

    ## fit the remaining nuisance paramters for time-to-event outcomes
    if (outcome_type == "time-to-event") {
      failure_hazard_fit <- fit_failure_hazard(
        train_data = data,
        valid_data = NULL,
        learners = failure_hazard_estimator,
        exposure = exposure,
        confounders = confounders,
        times = outcome
      )
      censoring_hazard_fit <- fit_censoring_hazard(
        train_data = data,
        valid_data = NULL,
        learners = censoring_hazard_estimator,
        exposure = exposure,
        confounders = confounders,
        times = outcome,
        censoring = censoring
      )
      cond_outcome_fit <- NULL
    }

    ## compute the efficient influence function
    ueif_dt <- uncentered_eif(
      data = data,
      scale = param_scale,
      confounders = confounders,
      exposure = exposure,
      outcome = outcome,
      modifiers = modifiers,
      prop_score_fit = prop_score_fit,
      prop_score_values = prop_score_values,
      cond_outcome_fit = cond_outcome_fit,
      failure_hazard_fit = failure_hazard_fit,
      censoring_hazard_fit = censoring_hazard_fit
    )

    ## estimate the estimands
    if (estimator == "onestep") {
      tem_vip_fit <- one_step_estimator(uncentered_eif_data = ueif_dt)
    } else {
      tem_vip_fit <- tml_estimator(
        data = data,
        confounders = confounders,
        modifiers = modifiers,
        exposure = exposure,
        outcome = outcome,
        scale = param_scale,
        prop_score_fit = prop_score_fit,
        cond_outcome_fit = cond_outcome_fit,
        failure_hazard_fit = failure_hazard_fit,
        censoring_hazard_fit = censoring_hazard_fit
      )
    }
  } else {

    ## make the validation folds
    if (outcome_type == "time-to-event") {
      folds <- origami::make_folds(
        data, V = cross_fit_folds, cluster_ids = data$id
      )
    } else if (outcome_type %in% c("continuous", "binary")) {
      folds <- origami::make_folds(data, V = cross_fit_folds)
    }

    ## perform cross-fitting procedure
    results <- origami::cross_validate(
      cv_fun = cross_fit_fold,
      folds = folds,
      data = data,
      confounders = confounders,
      modifiers = modifiers,
      exposure = exposure,
      outcome = outcome,
      censoring = censoring,
      outcome_type = outcome_type,
      scale = param_scale,
      estimator = estimator,
      cond_outcome_estimator = cond_outcome_estimator,
      prop_score_estimator = prop_score_estimator,
      prop_score_values = prop_score_values,
      failure_hazard_estimator = failure_hazard_estimator,
      censoring_hazard_estimator = censoring_hazard_estimator,
      .combine = FALSE,
      use_future = parallel
    )

    ## combine the estimates
    weighted_estimates <- lapply(
      seq_len(cross_fit_folds),
      function(idx) results$tem_vip_fit[[idx]] * results$prop_valid_data[[idx]]
    )
    weighted_estimates_dt <- data.table::rbindlist(weighted_estimates)
    tem_vip_fit <- weighted_estimates_dt[, lapply(.SD, sum)]
  }

  ## compute the confidence intervals and rescale everything
  if (cross_fit) {
    unscaled_eif_vars_dt <- data.table::rbindlist(results$unscaled_eif_vars)
    eif_vars <- unscaled_eif_vars_dt[, lapply(.SD, sum)] / n_obs
  } else {
    eif_vars <- ueif_dt[, lapply(.SD, var)]
  }
  test_dt <- test_hypotheses(
    n_obs = n_obs,
    estimates = tem_vip_fit,
    var_estimates = eif_vars,
    rescale_factor = rescale_factor
  )

  ## organize table in decreasing order of p-value
  test_dt <- test_dt[order(p_value), ]

  return(test_dt)
}



#' @title Cross-Fitting Procedure
#'
#' @description \code{cross_fit_fold()} estimates the treatment effect
#'   modification variable importance parameters using cross-fitting.
#'
#' @param fold An \code{\link[origami:make_folds]{origami}} fold object.
#' @param data A \code{data.table} containing the observed data.
#' @param confounders A \code{character} vector of column names corresponding to
#'   baseline covariates.
#' @param modifiers A \code{character} vector of columns names corresponding to
#'   the suspected effect modifiers. This vector must be a subset of
#'   \code{confounders}.
#' @param exposure A \code{character} corresponding to the exposure variable.
#' @param outcome A \code{character} corresponding to the outcome variable.
#' @param censoring A \code{character} indicating the right censoring indicator
#'   variable. Only used with time-to-event outcomes. Defaults to \code{NULL}.
#' @param outcome_type A \code{character} indicating the outcome type.
#'   \code{"continuous"}, \code{"binary"} and \code{"time-to-event"} are
#'   currently supported.
#' @param scale A \code{character} indicating the type of treatment effect
#'   modifier variable importance parameter. Currently supports
#'   \code{"absolute"} and \code{"relative"}.
#' @param estimator A \code{character} set to either \code{"tmle"} or
#'   \code{"onestep"}. The former results in \code{unihtee()} to use a targeted
#'   maximum likelihood estimators to estimate the deisred TEM VIP, while the
#'   latter uses a one step estimator.
#' @param cond_outcome_estimator A \code{\link[sl3]{Stack}}, or other learner
#'   class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
#'   learners from \pkg{sl3} to estimate the conditional outcome. Defaults to a
#'   generalized linear model with one- and two- way interactions among all
#'   \code{confounders} and \code{exposure} variables. Only used with continuous
#'   and binary outcomes.
#' @param prop_score_estimator A \code{\link[sl3]{Stack}}, or other learner
#'   class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
#'   learners from \pkg{sl3} to estimate the propensity score. Defaults to a
#'   generalized linear model with one- and two- way interactions among all
#'   \code{confounders} variables.
#' @param prop_score_values An optional \code{character} corresponding to the
#'   (known) propensity score values for each observation in \code{data}.
#'   Defaults to \code{NULL}.
#' @param failure_hazard_estimator A \code{\link[sl3]{Stack}}, or other learner
#'   class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a set of
#'   learners from \pkg{sl3} to estimate the conditional failure hazard
#'   function. Defaults to an xgboost learner with \code{confounders} and
#'   \code{exposure} variables as covariates. Only used with time-to-event
#'   outcomes.
#' @param censoring_hazard_estimator A \code{\link[sl3]{Stack}}, or other
#'   learner class (inheriting from \code{\link[sl3]{Lrnr_base}}), containing a
#'   set of learners from \pkg{sl3} to estimate the conditional censoring hazard
#'   function. Defaults to an xgboost learner with \code{confounders} and
#'   \code{exposure} variables as covariates. Only used with time-to-event
#'   outcomes.
#'
#' @return A \code{list} object containing the validation dataset's uncentered
#'   efficient influence function estimates, the treatment effect modification
#'   variable importance parameter estimates and the proportion of observations
#'   in the validation data.
#'
#' @importFrom data.table as.data.table .I .SD
#' @importFrom origami training validation
#'
#' @keywords internal
#'
cross_fit_fold <- function(fold,
                           data,
                           confounders,
                           modifiers,
                           exposure,
                           outcome,
                           censoring,
                           outcome_type,
                           scale,
                           estimator,
                           cond_outcome_estimator,
                           prop_score_estimator,
                           prop_score_values,
                           failure_hazard_estimator,
                           censoring_hazard_estimator) {

  ## split the data into training and testing
  train_data <- origami::training(data)
  valid_data <- origami::validation(data)

  ## fit the propensity score, if necessary
  if (outcome_type == "time-to-event") {
    prop_score_train_data <- train_data[train_data[, .I[1], by = "id"]$V1]
  } else {
    prop_score_train_data <- data.table::copy(train_data)
  }
  if (is.null(prop_score_values)) {
    prop_score_fit <- fit_prop_score(
      train_data = prop_score_train_data,
      valid_data = valid_data,
      learners = prop_score_estimator,
      exposure = exposure,
      confounders = confounders
    )
  } else {
    prop_score_fit <- list(
      estimates = data[[prop_score_values]]
    )
  }

  ## fit the remaining nuisance parameter for continuous and binary outcomes
  if (outcome_type %in% c("continuous", "binary")) {
    cond_outcome_fit <- fit_cond_outcome(
      train_data = train_data,
      valid_data = valid_data,
      learners = cond_outcome_estimator,
      exposure = exposure,
      confounders = confounders,
      outcome = outcome
    )
    failure_hazard_fit <- NULL
    censoring_hazard_fit <- NULL
  }

  ## fit the remaining nuisance paramters for time-to-event outcomes
  if (outcome_type == "time-to-event") {
    failure_hazard_fit <- fit_failure_hazard(
      train_data = train_data,
      valid_data = valid_data,
      learners = failure_hazard_estimator,
      exposure = exposure,
      confounders = confounders,
      times = outcome
    )
    censoring_hazard_fit <- fit_censoring_hazard(
      train_data = train_data,
      valid_data = valid_data,
      learners = censoring_hazard_estimator,
      exposure = exposure,
      confounders = confounders,
      times = outcome,
      censoring = censoring
    )
    cond_outcome_fit <- NULL
  }

  ## compute the efficient influence function
  ueif_dt <- uncentered_eif(
    data = valid_data,
    scale = scale,
    confounders = confounders,
    exposure = exposure,
    outcome = outcome,
    modifiers = modifiers,
    prop_score_fit = prop_score_fit,
    prop_score_values = prop_score_values,
    cond_outcome_fit = cond_outcome_fit,
    failure_hazard_fit = failure_hazard_fit,
    censoring_hazard_fit = censoring_hazard_fit
  )

  ## estimate the estimands
  if (estimator == "onestep") {
    tem_vip_fit <- one_step_estimator(uncentered_eif_data = ueif_dt)
  } else {
    tem_vip_fit <- tml_estimator(
      data = valid_data,
      confounders = confounders,
      modifiers = modifiers,
      exposure = exposure,
      outcome = outcome,
      scale = scale,
      prop_score_fit = prop_score_fit,
      cond_outcome_fit = cond_outcome_fit,
      failure_hazard_fit = failure_hazard_fit,
      censoring_hazard_fit = censoring_hazard_fit
    )
  }

  ## return the centered efficient influence function, estimates and
  ## proportion of observations in the validation data
  unscaled_eif_vars <- nrow(ueif_dt) * ueif_dt[, lapply(.SD, var)]
  return(list(
    "unscaled_eif_vars" = unscaled_eif_vars,
    "tem_vip_fit" = tem_vip_fit,
    "prop_valid_data" = nrow(valid_data) / nrow(data)
  ))
}
