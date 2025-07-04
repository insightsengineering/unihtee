source("testing-utils.R")

test_that(
  paste("one_step_estimator() produces accurate estimates without",
        "cross-fitting (cont outcome)"),
  {
    library(sl3)

    # generate data
    set.seed(72342)
    dt <- generate_test_data(n_obs = 1000)
    y_min <- min(dt$y)
    y_max <- max(dt$y)
    dt$y_scaled <- (dt$y - y_min) / (y_max - y_min)

    # fit the propensity score
    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected cond outcome
    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_xgboost$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y_scaled"
    )

    # compute the uncentered eif
    ueif_dt <- uncentered_eif(
      data = dt,
      effect = "absolute",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y_scaled",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = prop_score_fit,
      prop_score_values = NULL,
      cond_outcome_fit = cond_outcome_fit,
      failure_hazard_fit = NULL,
      censoring_hazard_fit = NULL
    )

    one_step_fit <- one_step_estimator(uncentered_eif_data = ueif_dt)

    # note that the true parameter values for w_1, w_3 are 0, 1
    one_step_fit <- one_step_fit * (y_max - y_min)
    expect_equal(one_step_fit$w_1, 0, tolerance = 0.1)
    expect_equal(one_step_fit$w_3, 1, tolerance = 0.1)
  }
)

test_that(
  "one_step_estimator() solves the efficient influence function (cont outcome)",
  {
    library(sl3)

    # generate data
    set.seed(84891)
    dt <- generate_test_data(n_obs = 100)
    y_min <- min(dt$y)
    y_max <- max(dt$y)
    dt$y_scaled <- (dt$y - y_min) / (y_max - y_min)

    # fit the propensity score
    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected cond outcome
    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_xgboost$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y_scaled"
    )

    # compute the uncentered eif
    ueif_dt <- uncentered_eif(
      data = dt,
      effect = "absolute",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y_scaled",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = prop_score_fit,
      prop_score_values = NULL,
      cond_outcome_fit = cond_outcome_fit,
      failure_hazard_fit = NULL,
      censoring_hazard_fit = NULL
    )

    one_step_fit <- one_step_estimator(uncentered_eif_data = ueif_dt)
    one_step_fit <- one_step_fit * (y_max - y_min)
    ueif_dt <- ueif_dt * (y_max - y_min)

    # note that the true parameter values for w_1, w_3 are 0, 1
    expect_equal(
      mean(ueif_dt$w_1 - one_step_fit$w_1), 0, tolerance = 1e-5
    )
    expect_equal(
      mean(ueif_dt$w_3 - one_step_fit$w_3), 0, tolerance = 1e-5
    )
  }
)

test_that(
  paste("tml_estimator() produces accurate estimates without cross-fitting",
        "(cont outcome)"),
  {
    library(sl3)

    # generate data
    set.seed(5123)
    dt <- generate_test_data(n_obs = 500)
    y_min <- min(dt$y)
    y_max <- max(dt$y)
    dt$y_scaled <- (dt$y - y_min) / (y_max - y_min)

    # fit the propensity score
    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected cond outcome
    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_xgboost$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y_scaled"
    )

    # compute the TML estimate
    tmle_fit <- tml_estimator(
      data = dt,
      effect = "absolute",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y_scaled",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = prop_score_fit,
      cond_outcome_fit = cond_outcome_fit
    )

    # note that the true parameter values for w_1, w_3 are 0, 1
    tmle_fit <- tmle_fit * (y_max - y_min)
    expect_equal(tmle_fit$w_1, 0, tolerance = 0.1)
    expect_equal(tmle_fit$w_3, 1, tolerance = 0.1)
  }
)

test_that(
  paste("one_step_estimator() solves the efficient influence function",
        "(time-to-event outcome, RD TEM VIP)"),
  {
    library(sl3)

    # generate data
    set.seed(84891)
    dt <- generate_test_data(n_obs = 100, outcome_type = "time-to-event")
    long_dt <- tte_data_melt(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = 5,
      prop_score_values = NULL
    )

    # fit the propensity score
    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = long_dt,
      learners = sl3::Lrnr_xgboost$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected failure hazard
    fail_fit <- fit_failure_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3:::Lrnr_xgboost$new(),
      exposure = "a",
      times = "time",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected failure hazard
    cens_fit <- fit_censoring_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3:::Lrnr_xgboost$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      times = "time",
      censoring = "censoring"
    )

    # compute the uncentered eif
    ueif_dt <- uncentered_eif(
      data = long_dt,
      effect = "absolute",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = prop_score_fit,
      cond_outcome_fit = NULL,
      prop_score_values = NULL,
      failure_hazard_fit = fail_fit,
      censoring_hazard_fit = cens_fit
    )


    one_step_fit <- one_step_estimator(uncentered_eif_data = ueif_dt)

    # note that the true parameter values for w_1, w_3 are 0, 1
    expect_equal(c(
      mean(ueif_dt$w_1 - one_step_fit$w_1),
      mean(ueif_dt$w_3 - one_step_fit$w_3)
    ),
    c(0, 0),
    tolerance = 1e-10
    )
  }
)

test_that(
  paste("one_step_estimator() produces accurate estimates without",
        "cross-fitting (time-to-event outcome, RD TEM VIP)"),
  {
    library(sl3)

    # generate data
    set.seed(84891)
    dt <- generate_test_data(n_obs = 5000, outcome_type = "time-to-event")
    long_dt <- tte_data_melt(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = 5,
      prop_score_values = NULL
    )

    # fit the propensity score
    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = long_dt,
      learners = sl3::Lrnr_xgboost$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected failure hazard
    fail_fit <- fit_failure_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3:::Lrnr_xgboost$new(),
      exposure = "a",
      times = "time",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected failure hazard
    cens_fit <- fit_censoring_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3:::Lrnr_xgboost$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      times = "time",
      censoring = "censoring"
    )

    # compute the uncentered eif
    ueif_dt <- uncentered_eif(
      data = long_dt,
      effect = "absolute",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      modifiers = c("w_1", "w_2", "w_3"),
      prop_score_fit = prop_score_fit,
      cond_outcome_fit = NULL,
      prop_score_values = NULL,
      failure_hazard_fit = fail_fit,
      censoring_hazard_fit = cens_fit
    )

    # fit the propensity score
    one_step_fit <- one_step_estimator(uncentered_eif_data = ueif_dt)

    # note that the true parameter values for w_1, w_3 are approx 1.7, 0, 0
    expect_equal(one_step_fit$w_1, 1.7, tolerance = 0.1)
    expect_equal(one_step_fit$w_2, 0, tolerance = 0.1)
    expect_equal(one_step_fit$w_3, 0, tolerance = 0.1)
  }
)

test_that(
  paste("tml_estimator() produces accurate estimates without",
        "cross-fitting (time-to-event outcome, RD TEM VIP)"),
  {
    library(sl3)

    # generate data
    set.seed(84891)
    dt <- generate_test_data(n_obs = 5000, outcome_type = "time-to-event")
    long_dt <- tte_data_melt(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = 5,
      prop_score_values = NULL
    )

    # fit the propensity score
    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = long_dt,
      learners = sl3::Lrnr_xgboost$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected failure hazard
    fail_fit <- fit_failure_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3::Lrnr_xgboost$new(),
      exposure = "a",
      times = "time",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected failure hazard
    cens_fit <- fit_censoring_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3:::Lrnr_glmnet$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      times = "time",
      censoring = "censoring"
    )

    # compute the TML estimate
    tmle_fit <- tml_estimator(
      data = long_dt,
      confounders = c("w_1", "w_2", "w_3"),
      modifiers = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      effect = "absolute",
      prop_score_values = NULL,
      prop_score_fit = prop_score_fit,
      cond_outcome_fit = NULL,
      failure_hazard_fit = fail_fit,
      censoring_hazard_fit = cens_fit
    )

    # note that the true parameter values for w_1, w_3 are approx 1.7, 0, 0
    expect_equal(tmle_fit$w_1, 1.72, tolerance = 0.1)
    expect_equal(tmle_fit$w_2, 0, tolerance = 0.1)
    expect_equal(tmle_fit$w_3, 0, tolerance = 0.1)
  }
)


test_that(
  paste("plugin_estimator() for continuous outcomes don't have enormous error"),
  {
    library(sl3)

    # generate data
    set.seed(72342)
    dt <- generate_test_data(n_obs = 100000)

    # fit the expected cond outcome
    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_xgboost$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    # plugin estimates
    plugin_estimates <- plugin_estimator(
      data = dt,
      effect = "absolute",
      outcome = "y",
      modifiers = c("w_1", "w_3"),
      cond_outcome_fit = cond_outcome_fit,
      failure_hazard_fit = NULL,
      censoring_hazard_fit = NULL
    )

    # make sure error isn't too egregious asymptotically
    expect_equal(as.numeric(plugin_estimates["w_1"]), 0, tolerance = 0.1)
    expect_equal(as.numeric(plugin_estimates["w_3"]) - 1, 0, tolerance = 0.1)
  }
)

test_that(
  paste("plugin_estimator() for time-to-event outcomes don't have enormous",
        "error"),
  {
    library(sl3)

    # generate data
    set.seed(72342)
    dt <- generate_test_data(n_obs = 5000, outcome_type = "time-to-event")
    long_dt <- tte_data_melt(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = 5,
      prop_score_values = NULL
    )

    # fit the expected failure hazard
    fail_fit <- fit_failure_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3::Lrnr_xgboost$new(),
      exposure = "a",
      times = "time",
      confounders = c("w_1", "w_2", "w_3")
    )

    # plugin estimates
    plugin_estimates <- plugin_estimator(
      data = long_dt,
      effect = "absolute",
      outcome = "time",
      modifiers = c("w_1", "w_2", "w_3"),
      cond_outcome_fit = NULL,
      failure_hazard_fit = fail_fit
    )

    # make sure error isn't too egregious asymptotically
    expect_equal(as.numeric(plugin_estimates["w_1"]) - 1.72, 0, tolerance = 0.2)
    expect_equal(as.numeric(plugin_estimates["w_2"]), 0, tolerance = 0.2)
    expect_equal(as.numeric(plugin_estimates["w_3"]), 0, tolerance = 0.2)
  }
)
