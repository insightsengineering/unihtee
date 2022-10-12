source("testing-utils.R")

test_that(
  "one_step_estimator() produces accurate estimates without cross-fitting", {

    library(sl3)

    # generate data
    set.seed(84891)
    dt <- generate_test_data(n_obs = 5000)

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
      learners = sl3::Lrnr_ranger$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    # compute the uncentered eif
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

    one_step_fit <- one_step_estimator(uncentered_eif_data = ueif_dt)

    # note that the true parameter values for w_1, w_3 are 0, 1
    expect_equal(as.numeric(one_step_fit), c(0, 1), tolerance = 0.1)
})

test_that(
  "one_step_estimator() solves the efficient influence function", {

    library(sl3)

    # generate data
    set.seed(84891)
    dt <- generate_test_data(n_obs = 100)

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
      learners = sl3::Lrnr_ranger$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    # compute the uncentered eif
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

    one_step_fit <- one_step_estimator(uncentered_eif_data = ueif_dt)

    # note that the true parameter values for w_1, w_3 are 0, 1
    expect_equal(c(mean(ueif_dt$w_1 - one_step_fit$w_1),
                   mean(ueif_dt$w_3 - one_step_fit$w_3)),
                 c(0, 0),
                 tolerance = 1e-10)
})
