source("testing-utils.R")

test_that(
  paste(
    "uncentered_eif() minus true risk diff parameter value has a mean",
    "of zero when propensity scores aren't known"
  ),
  {
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
      type = "risk difference",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = prop_score_fit,
      prop_score_values = NULL,
      cond_outcome_fit = cond_outcome_fit
    )

    # note that the true parameter values are equal to 0,1 for W_1, W_3
    expect_equal(eif[, sapply(.SD, mean)], c("w_1" = 0, "w_3" = 1),
      tolerance = 0.005
    )
  }
)

test_that(
  paste(
    "uncentered_eif() minus true risk diff parameter value has a mean",
    "of zero when propensity scores are known"
  ),
  {
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
      type = "risk difference",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = NULL,
      cond_outcome_fit = cond_outcome_fit,
      prop_score_values = dt$prop_score
    )

    # note that the true parameter values are equal to 0,1 for W_1, W_3
    expect_equal(eif[, sapply(.SD, mean)], c("w_1" = 0, "w_3" = 1),
      tolerance = 0.01
    )
  }
)

test_that(
  paste(
    "uncentered_eif() minus true risk diff parameter value has a mean",
    "of zero when propensity scores aren't known"
  ),
  {
    library(sl3)

    # generate data
    set.seed(5123)
    dt <- generate_test_data(n_obs = 100000, outcome_type = "binary")

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
      type = "relative risk",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = prop_score_fit,
      prop_score_values = NULL,
      cond_outcome_fit = cond_outcome_fit
    )

    # note that the true parameter values are approx equal to 0, 4.22 for W_1,
    # W_3
    expect_equal(eif[, sapply(.SD, mean)], c("w_1" = 0, "w_3" = 4.2),
      tolerance = 0.1
    )
  }
)

test_that(
  paste(
    "uncentered_eif() minus true risk diff parameter value has a mean",
    "of zero when propensity scores are known"
  ),
  {
    library(sl3)

    # generate data
    set.seed(5123)
    dt <- generate_test_data(n_obs = 100000, outcome_type = "binary")

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
      type = "relative risk",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = NULL,
      cond_outcome_fit = cond_outcome_fit,
      prop_score_values = dt$prop_score
    )

    # note that the true parameter values are equal approx to 0, 4.22 for W_1,
    # W_3
    expect_equal(eif[, sapply(.SD, mean)], c("w_1" = 0, "w_3" = 4.2),
      tolerance = 0.1
    )
  }
)
