source("testing-utils.R")

test_that(
  "one_step_estimator() produces accurate estimates without cross-fitting",
  {
    library(sl3)

    # generate data
    set.seed(234)
    dt <- generate_test_data(n_obs = 1000, outcome_type = "binary")

    # fit the propensity score
    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm$new(),
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
      outcome = "y"
    )

    # compute the uncentered eif
    ueif_dt <- uncentered_eif(
      data = dt,
      effect = "relative",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = prop_score_fit,
      prop_score_values = NULL,
      cond_outcome_fit = cond_outcome_fit,
      failure_hazard_fit = NULL,
      censoring_hazard_fit = NULL
    )

    one_step_fit <- one_step_estimator(uncentered_eif_data = ueif_dt)

    # note that the true parameter values for w_3 is 2.0
    expect_equal(one_step_fit$w_3, 2.0, tolerance = 0.1)
  }
)

test_that(
  "one_step_estimator() solves the efficient influence function",
  {
    library(sl3)

    # generate data
    set.seed(84891)
    dt <- generate_test_data(n_obs = 100, outcome_type = "binary")

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
      learners = sl3::Lrnr_glmnet$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    # compute the uncentered eif
    ueif_dt <- uncentered_eif(
      data = dt,
      effect = "relative",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = prop_score_fit,
      prop_score_values = NULL,
      cond_outcome_fit = cond_outcome_fit,
      failure_hazard_fit = NULL,
      censoring_hazard_fit = NULL
    )

    one_step_fit <- one_step_estimator(uncentered_eif_data = ueif_dt)

    expect_equal(mean(ueif_dt$w_1 - one_step_fit$w_1), 0, tolerance = 1e-5)
    expect_equal(mean(ueif_dt$w_3 - one_step_fit$w_3), 0, tolerance = 1e-5)
  }
)

test_that(
  "tml_estimator() produces accurate estimates without cross-fitting",
  {
    library(sl3)

    # generate data
    set.seed(514)
    dt <- generate_test_data(n_obs = 5000, outcome_type = "binary")

    # fit the propensity score
    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glmnet$new(family = "binomial"),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected cond outcome
    lrnr_lasso <- sl3::Lrnr_glmnet$new(
      family = "binomial",
      formula = "~ a * w_1 + a * w_2 + a * w_3"
    )
    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = lrnr_lasso,
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    # compute the TML estimate
    tmle_fit <- tml_estimator(
      data = dt,
      effect = "relative",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = prop_score_fit,
      cond_outcome_fit = cond_outcome_fit
    )

    # note that the true parameter values for w_1, w_3 are 0, 2.0
    expect_equal(tmle_fit$w_1, 0, tolerance = 0.1)
    expect_equal(tmle_fit$w_3, 2.0, tolerance = 0.1)
  }
)

test_that(
  paste("one_step_estimator() produces accurate estimates without",
        "cross-fitting, TTE outcome"),
  {
    # generate data
    set.seed(426723)
    dt <- generate_test_data(n_obs = 5000, outcome_type = "time-to-event RR")
    long_dt <- tte_data_melt(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = 3,
      prop_score_values = NULL
    )

    # fit the propensity score
    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = long_dt,
      learners = sl3::Lrnr_glm$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected failure hazard
    fail_fit <- fit_failure_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3:::Lrnr_earth$new(
        formula = "~ a * w_1 + w_2 + w_3",
        glm = list(family = "binomial")
      ),
      exposure = "a",
      times = "time",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected failure hazard
    cens_fit <- fit_censoring_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3:::Lrnr_glm$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      times = "time",
      censoring = "censoring"
    )

    eif <- uncentered_eif(
      data = long_dt,
      effect = "relative",
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

    # fit the estimator
    one_step_fit <- one_step_estimator(uncentered_eif_data = eif)

    # note that the true parameter values for w_1, w_2, w_3 are 0.6, 0, 0
    expect_equal(abs(one_step_fit$w_1 - 0.6), 0, tolerance = 0.1)
    expect_equal(one_step_fit$w_2, 0, tolerance = 0.1)
    expect_equal(one_step_fit$w_3, 0, tolerance = 0.1)
  }
)

test_that(
  paste("one_step_estimator() solves the efficient influence function, TTE",
        "outcome"),
  {
    library(sl3)

    # generate data
    set.seed(42)
    dt <- generate_test_data(n_obs = 100, outcome_type = "time-to-event RR")
    long_dt <- tte_data_melt(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = 3,
      prop_score_values = "prop_score"
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
      learners = sl3:::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      times = "time",
      censoring = "censoring"
    )

    eif <- uncentered_eif(
      data = long_dt,
      effect = "relative",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      modifiers = c("w_1", "w_2", "w_3"),
      prop_score_fit = NULL,
      cond_outcome_fit = NULL,
      prop_score_values = "prop_score",
      failure_hazard_fit = fail_fit,
      censoring_hazard_fit = cens_fit
    )

    # fit the estimator
    one_step_fit <- one_step_estimator(uncentered_eif_data = eif)

    expect_equal(c(
      mean(eif$w_1 - one_step_fit$w_1),
      mean(eif$w_2 - one_step_fit$w_2),
      mean(eif$w_3 - one_step_fit$w_3)
    ),
    c(0, 0, 0),
    tolerance = 1e-10
    )
  }
)

test_that(
  paste("tml_estimator() produces accurate estimates without",
        "cross-fitting, TTE outcome"),
  {

    library(earth)

    # generate data
    set.seed(72342)
    dt <- generate_test_data(n_obs = 2000, outcome_type = "time-to-event RR")
    long_dt <- tte_data_melt(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = 3,
      prop_score_values = NULL
    )

    # fit the propensity score
    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = long_dt,
      learners = sl3::Lrnr_glm$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected failure hazard
    fail_fit <- fit_failure_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3:::Lrnr_earth$new(
        formula = "~ a * w_1 + w_2 + w_3",
        glm = list(family = "binomial")
      ),
      exposure = "a",
      times = "time",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected failure hazard
    cens_fit <- fit_censoring_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3:::Lrnr_glm$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      times = "time",
      censoring = "censoring"
    )

    # TML estimates
    tmle <- tml_estimator(
      data = long_dt,
      confounders = c("w_1", "w_2", "w_3"),
      modifiers = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      effect = "relative",
      prop_score_values = NULL,
      prop_score_fit = prop_score_fit,
      cond_outcome_fit = NULL,
      failure_hazard_fit = fail_fit,
      censoring_hazard_fit = cens_fit
    )

    # note that the true parameter values for w_1, w_2, w_3 are 0.6, 0, 0
    expect_equal(abs(tmle$w_1 - 0.6), 0, tolerance = 0.1)
    expect_equal(tmle$w_2, 0, tolerance = 0.1)
    expect_equal(tmle$w_3, 0, tolerance = 0.1)
  }
)
