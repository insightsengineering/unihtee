source("testing-utils.R")

test_that(
  paste(
    "fit_prop_score() outputs a vector of estimates that is similar to",
    "ground truth when a validation dataset isn't provided"
  ),
  {
    library(sl3)

    # generate data
    set.seed(5123)
    dt <- generate_test_data(n_obs = 200)

    # fit the propensity score
    fit <- fit_prop_score(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(
        formula = "~ w_1 + w_2 + w_3"
      ),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    # make sure that the MSE is approximately zero
    expect_equal(mean((fit$estimates - dt$prop_score)^2), 0, tolerance = 0.01)
  }
)

test_that(
  "fit_prop_score() outputs the sl3 lrnr object fit on the training data",
  {
    library(sl3)

    # generate data
    set.seed(5123)
    dt <- generate_test_data(n_obs = 200)

    # fit the propensity score
    fit <- fit_prop_score(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(
        formula = "~ w_1 + w_2 + w_3"
      ),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    expect_equal(class(fit$fit), c("Lrnr_glm_fast", "Lrnr_base", "R6"))
  }
)

test_that(
  paste(
    "fit_prop_score() outputs a vector of estimates that is similar to",
    "ground truth when a validation set is provided"
  ),
  {
    library(sl3)

    # generate data
    set.seed(5123)
    train_dt <- generate_test_data(n_obs = 400)
    valid_dt <- generate_test_data(n_obs = 200)

    # fit the propensity score
    interactions <- list(c("w_1", "a"))
    lrnr_enet <- sl3::Lrnr_glmnet$new(
      alpha = 0.5,
      formula = "~w_1 + w_2 + w_3"
    )
    fit <- fit_prop_score(
      train_data = train_dt,
      valid_data = valid_dt,
      learners = lrnr_enet,
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    # make sure that the MSE is approximately zero
    expect_equal(mean((fit$estimates - valid_dt$prop_score)^2), 0,
      tolerance = 0.01
    )
  }
)

test_that(
  paste(
    "fit_cond_est() returns a vector of estimates that is similar to",
    "ground truth when no validation set is provided"
  ),
  {
    library(sl3)

    # generate data
    set.seed(12312)
    dt <- generate_test_data(n_obs = 1000)

    # fit the expected cond outcome
    fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_ranger$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    # make sure that the MSE is approximately zero
    expect_equal(mean((fit$estimates - dt$y)^2), 0, tolerance = 0.05)
  }
)

test_that(
  paste(
    "fit_cond_est() returns a vector of estimates that is similar to",
    "ground truth when a validation set is provided"
  ),
  {
    library(sl3)

    # generate data
    set.seed(62412)
    train_dt <- generate_test_data(n_obs = 5000)
    valid_dt <- generate_test_data(n_obs = 200)

    # fit the expected cond outcome
    fit <- fit_cond_outcome(
      train_data = train_dt,
      valid_data = valid_dt,
      learners = sl3::Lrnr_ranger$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    # make sure that the MSE is approximately zero
    expect_equal(mean((fit$estimates - valid_dt$y)^2), 0, tolerance = 0.05)
  }
)

test_that(
  "fit_cond_outcome() outputs the sl3 lrnr object fit on the training data",
  {
    library(sl3)

    # generate data
    set.seed(5123)
    dt <- generate_test_data(n_obs = 200)

    # fit the expected cond outcome
    fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    expect_equal(class(fit$fit), c("Lrnr_glm_fast", "Lrnr_base", "R6"))
  }
)

test_that(
  paste(
    "fit_cond_est() returns vectors of potential outcome estimates that is",
    "similar to ground truth when no validation set is provided"
  ),
  {
    library(sl3)

    # generate data
    set.seed(12312)
    dt <- generate_test_data(n_obs = 5000)

    # fit the expected cond outcome
    fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_ranger$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    # make sure that the MSE is approximately zero
    expect_equal(mean((fit$exp_estimates - dt$y_1)^2), 0, tolerance = 0.1)
    expect_equal(mean((fit$noexp_estimates - dt$y_0)^2), 0, tolerance = 0.1)
  }
)

test_that(
  paste(
    "fit_failure_hazard() returns vectors of potential conditional failure",
    "hazard estmiates that are close to the ground truth"
  ),
  {
    library(sl3)

    # generate the data
    set.seed(100)
    dt <- generate_test_data(n_obs = 300, outcome_type = "time-to-event")
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
    fit <- fit_failure_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3:::Lrnr_ranger$new(),
      exposure = "a",
      times = "time",
      confounders = c("w_1", "w_2", "w_3")
    )

    # compute the true hazards
    cond_surv_hazard <- function(time, exposure, w_1, w_2, w_3) {
      (time < 9) / (1 + exp(2 + 3 * exposure * w_1)) + (time == 9)
    }
    exp_truth <-
      cond_surv_hazard(long_dt$time, 1, long_dt$w_1, long_dt$w_2, long_dt$w_3)
    noexp_truth <-
      cond_surv_hazard(long_dt$time, 0, long_dt$w_1, long_dt$w_2, long_dt$w_3)

    expect_equal(mean((fit$exp_estimates - exp_truth)^2), 0, tolerance = 0.05)
    expect_equal(
      mean((fit$noexp_estimates - noexp_truth)^2), 0, tolerance = 0.05
    )
  }
)

test_that(
  paste(
    "fit_censoring_hazard() returns vectors of potential conditional censoring",
    "hazard estimates that are close to the ground truth"
  ),
  {
    library(sl3)

    # generate the data
    set.seed(100)
    dt <- generate_test_data(n_obs = 300, outcome_type = "time-to-event")
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
    fit <- fit_censoring_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3:::Lrnr_glm_fast$new(
        formula = "~ w_1 * a + w_2 * a + w_3 * a"
      ),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      times = "time",
      censoring = "censoring"
    )

    # compute the true hazards
    expect_equal(mean((fit$exp_estimates - 0.05)^2), 0, tolerance = 0.005)
    expect_equal(
      mean((fit$noexp_estimates - 0.05)^2), 0, tolerance = 0.005
    )
  }
)

test_that("one_step_ate_estimator() estimates are close to gound truth", {

  # generate data
  set.seed(723423)
  dt <- generate_test_data(n_obs = 1000)

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
    learners = sl3::make_learner(
      sl3::Pipeline,
      sl3::Lrnr_define_interactions$new(
        list(c("w_1", "a"), c("w_2", "a"), c("w_3", "a"))
      ),
      sl3::Lrnr_glm_fast$new()
    ),
    exposure = "a",
    confounders = c("w_1", "w_2", "w_3"),
    outcome = "y"
  )

  # estimate the ATE
  one_step_ate_estimate <- one_step_ate_estimator(
    data = dt,
    confounders = c("w_1", "w_2", "w_3"),
    exposure = "a",
    outcome = "y",
    prop_score_fit = prop_score_fit,
    prop_score_values = NULL,
    cond_outcome_fit = cond_outcome_fit
  )

  # make sure that the ATE is near the true value of 1
  expect_equal(abs(one_step_ate_estimate - 1), 0, tolerance = 0.01)

})


test_that("tml_ate_estimator() estimates are close to gound truth", {

  # generate data
  set.seed(93453)
  dt <- generate_test_data(n_obs = 1000)

  # rescale the outomce
  outcome <- "y"
  min_out <- min(dt[[outcome]])
  max_out <- max(dt[[outcome]])
  rescale_factor <- max_out - min_out
  dt[[outcome]] <- (dt[[outcome]] - min_out) / rescale_factor

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
    learners = sl3::make_learner(
      sl3::Pipeline,
      sl3::Lrnr_define_interactions$new(
        list(c("w_1", "a"), c("w_2", "a"), c("w_3", "a"))
      ),
      sl3::Lrnr_glm_fast$new()
    ),
    exposure = "a",
    confounders = c("w_1", "w_2", "w_3"),
    outcome = "y"
  )

  # estimate the ATE
  tml_ate_estimate <- tml_ate_estimator(
    data = dt,
    confounders = c("w_1", "w_2", "w_3"),
    exposure = "a",
    outcome = "y",
    prop_score_fit = prop_score_fit,
    prop_score_values = NULL,
    cond_outcome_fit = cond_outcome_fit
  ) * rescale_factor

  # make sure that the ATE is near the true value of 1
  expect_equal(abs(tml_ate_estimate - 1), 0, tolerance = 0.01)

})
