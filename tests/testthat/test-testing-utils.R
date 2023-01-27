source("testing-utils.R")

test_that(
  "test_hypotheses() outputs a data.table with one row per modifier",
  {
    library(sl3)

    # generate data
    set.seed(6127)
    dt <- generate_test_data(n_obs = 200)

    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_ranger$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    ueif_dt <- uncentered_eif(
      data = dt,
      effect = "absolute",
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

    eif_vars <- ueif_dt[, lapply(.SD, var)]

    test_dt <- test_hypotheses(
      n_obs = nrow(dt),
      estimates = one_step_estimator(ueif_dt),
      var_estimates = eif_vars,
      rescale_factor = 1
    )

    expect_equal(nrow(test_dt), 2)
  }
)

test_that(
  paste(
    "test_hypotheses() contains modifier, estimate, se, z, p_value,",
    "ci_lower", "ci_upper", "p_value_fdr columns"
  ),
  {
    library(sl3)

    # generate data
    set.seed(6127)
    dt <- generate_test_data(n_obs = 200)

    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_ranger$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    ueif_dt <- uncentered_eif(
      data = dt,
      effect = "absolute",
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

    eif_vars <-
      ueif_dt[, lapply(.SD, var)]

    test_dt <- test_hypotheses(
      n_obs = nrow(dt),
      estimates = one_step_estimator(ueif_dt),
      var_estimates = eif_vars,
      rescale_factor = 1
    )

    expect_equal(
      colnames(test_dt),
      c(
        "modifier", "estimate", "se", "z", "p_value",
        "ci_lower", "ci_upper", "p_value_fdr"
      )
    )
  }
)

test_that(
  "test_hypotheses() properly rescales the estimates",
  {
    library(sl3)

    # generate data
    set.seed(62341)
    dt <- generate_test_data(n_obs = 10000)

    # rescale the outcome to be between 0 and 1
    min_y <- min(dt$y)
    max_y <- max(dt$y)
    dt$y <- (dt$y - min_y) / (max_y - min_y)

    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_ranger$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    ueif_dt <- uncentered_eif(
      data = dt,
      effect = "absolute",
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

    eif_vars <- ueif_dt[, lapply(.SD, var)]

    test_dt <- test_hypotheses(
      n_obs = nrow(dt),
      estimates = one_step_estimator(ueif_dt),
      var_estimates = eif_vars,
      rescale_factor = max_y - min_y
    )

    expect_equal(test_dt$estimate[1], 0, tolerance = 0.05)
    expect_equal(test_dt$estimate[2], 1, tolerance = 0.05)
  }
)

test_that(
  "test_hypotheses() properly rescales the variance estimates",
  {
    library(sl3)

    # generate data
    set.seed(62341)
    dt <- generate_test_data(n_obs = 1000)

    # rescale the outcome to be between 0 and 1
    min_y <- min(dt$y)
    max_y <- max(dt$y)
    dt$y <- (dt$y - min_y) / (max_y - min_y)

    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_ranger$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    ueif_dt <- uncentered_eif(
      data = dt,
      effect = "absolute",
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

    eif_vars <- ueif_dt[, lapply(.SD, var)]

    rescaled_test_dt <- test_hypotheses(
      n_obs = nrow(dt),
      estimates = one_step_estimator(ueif_dt),
      var_estimates = eif_vars,
      rescale_factor = max_y - min_y
    )

    # generate data
    set.seed(62341)
    dt <- generate_test_data(n_obs = 1000)

    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_ranger$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    ueif_dt <- uncentered_eif(
      data = dt,
      effect = "absolute",
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

    eif_vars <- ueif_dt[, lapply(.SD, var)]

    test_dt <- test_hypotheses(
      n_obs = nrow(dt),
      estimates = one_step_estimator(ueif_dt),
      var_estimates = eif_vars,
      rescale_factor = 1
    )

    expect_equal(rescaled_test_dt$se[1], test_dt$se[1], tolerance = 0.001)
    expect_equal(rescaled_test_dt$se[2], test_dt$se[2], tolerance = 0.001)
  }
)
