library(here)
source(here("tests", "testthat", "testing-utils.R"))

test_that(
  paste("uncentered_eif() minus true risk diff parameter value has a mean",
        "of zero when propensity scores aren't known"), {

    library(sl3)

    # generate data
    set.seed(5123)
    dt <- generate_test_data(n_obs = 100000)

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
      learners = sl3::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    # compute the uncentered eif
    eif <- uncentered_eif(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = "w_3",
      prop_score_fit = prop_score_fit,
      prop_score_values = NULL,
      cond_outcome_fit = cond_outcome_fit
    )

    # note that the true parameter value is equal to 1
    param <- 1

  expect_equal(mean(eif$w_3 - param), 0, tolerance = 0.005)
})

test_that(
  paste("uncentered_eif() minus true risk diff parameter value has a mean",
        "of zero when propensity scores are known"), {

    library(sl3)

    # generate data
    set.seed(5123)
    dt <- generate_test_data(n_obs = 100000)

    # fit the expected cond outcome
    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    # compute the uncentered eif
    eif <- uncentered_eif(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = "w_3",
      prop_score_fit = NULL,
      cond_outcome_fit = cond_outcome_fit,
      prop_score_values = dt$prop_score
    )

    # note that the true parameter value is equal to 1
    param <- 1

  expect_equal(mean(eif$w_3 - param), 0, tolerance = 0.005)
})
