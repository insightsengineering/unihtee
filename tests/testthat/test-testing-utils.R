library(here)
source(here("tests", "testthat", "testing-utils.R"))

test_that(
  "test_hypotheses() outputs a data.table with one row per modifier", {

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
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = prop_score_fit,
      prop_score_values = NULL,
      cond_outcome_fit = cond_outcome_fit
    )

    eif_vars <- ueif_dt[, lapply(.SD, var)]

    test_dt <- test_hypotheses(
      n_obs = nrow(dt),
      estimates = one_step_estimator(ueif_dt),
      var_estimates = eif_vars
    )

    expect_equal(nrow(test_dt), 2)
})

test_that(
  paste("test_hypotheses() contains modifier, estimate, se, z, p_value,",
        "ci_lower", "ci_upper", "p_value_fdr columns"), {

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
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = prop_score_fit,
      prop_score_values = NULL,
      cond_outcome_fit = cond_outcome_fit
    )

    eif_vars <- ueif_dt[, lapply(.SD, var)]

    test_dt <- test_hypotheses(
      n_obs = nrow(dt),
      estimates = one_step_estimator(ueif_dt),
      var_estimates = eif_vars
    )

    expect_equal(colnames(test_dt),
                 c("modifier", "estimate", "se", "z", "p_value",
                   "ci_lower", "ci_upper", "p_value_fdr"))
})
