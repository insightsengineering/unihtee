source("testing-utils.R")

test_that(
  paste("unihtee() uncovers treatment effect modifiers on the risk diff scale",
        "when outcomes are continuous with the one-step estimator"), {

  # generate data
  set.seed(712435)
  dt <- generate_test_data(n_obs = 5000)

  # apply unihtee
  results <- unihtee(
    data = dt,
    confounders = c("w_1", "w_2", "w_3"),
    modifiers = c("w_1", "w_2", "w_3"),
    exposure = "a",
    outcome = "y",
    outcome_type = "continuous",
    risk_type = "risk difference",
    estimator = "onestep"
  )

  # ensure that the adjusted p-value of w_3 is less than 0.05, and those of
  # w_1 and w_2 are above 0.05
  expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))

  # check that the estimates are reasonably close to the ground truth
  expect_equal(results$estimate[1], 1, tolerance = 0.1)
  expect_equal(results$estimate[2], 0, tolerance = 0.1)
  expect_equal(results$estimate[3], 0, tolerance = 0.1)
})

test_that(
  paste("unihtee() uncovers treatment effect modifiers on the risk diff scale",
        "when outcomes are continuous with the TML estimator"), {

  # generate data
  set.seed(712535)
  dt <- generate_test_data(n_obs = 5000)

  # apply unihtee
  results <- unihtee(
    data = dt,
    confounders = c("w_1", "w_2", "w_3"),
    modifiers = c("w_1", "w_2", "w_3"),
    exposure = "a",
    outcome = "y",
    outcome_type = "continuous",
    risk_type = "risk difference",
    estimator = "tmle",
    cond_outcome_estimator = sl3::Lrnr_ranger$new()
  )

  # ensure that the adjusted p-value of w_3 is less than 0.05, and those of
  # w_1 and w_2 are above 0.05
  expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))

  # check that the estimates are reasonably close to the ground truth
  expect_equal(results$estimate[1], 1, tolerance = 0.1)
  expect_equal(results$estimate[2], 0, tolerance = 0.1)
  expect_equal(results$estimate[3], 0, tolerance = 0.1)
})

test_that(
  paste("unihtee() uncovers treatment effect modifiers on the risk diff scale",
        "when outcomes are binary with the one-step estimator"), {

  # generate data
  set.seed(712435)
  dt <- generate_test_data(n_obs = 5000, outcome_type = "binary")

  # apply unihtee
  results <- unihtee(
    data = dt,
    confounders = c("w_1", "w_2", "w_3"),
    modifiers = c("w_1", "w_2", "w_3"),
    exposure = "a",
    outcome = "y",
    outcome_type = "binary",
    risk_type = "risk difference",
    estimator = "onestep"
  )

  # ensure that the adjusted p-value of w_3 is less than 0.05, and those of
  # w_1 and w_2 are above 0.05
  expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))

  # check that the estimates are reasonably close to the ground truth
  expect_equal(results$estimate[1], 0.38, tolerance = 0.1)
  expect_equal(results$estimate[2], 0, tolerance = 0.1)
  expect_equal(results$estimate[3], 0, tolerance = 0.1)
})

test_that(
  paste("unihtee() uncovers treatment effect modifiers on the risk diff scale",
        "when outcomes are binary with the TML estimator"), {

  # generate data
  set.seed(712535)
  dt <- generate_test_data(n_obs = 5000, outcome_type = "binary")

  # apply unihtee
  results <- unihtee(
    data = dt,
    confounders = c("w_1", "w_2", "w_3"),
    modifiers = c("w_1", "w_2", "w_3"),
    exposure = "a",
    outcome = "y",
    outcome_type = "binary",
    risk_type = "risk difference",
    estimator = "tmle"
  )

  # ensure that the adjusted p-value of w_3 is less than 0.05, and those of
  # w_1 and w_2 are above 0.05
  expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))

  # check that the estimates are reasonably close to the ground truth
  expect_equal(results$estimate[1], 0.38, tolerance = 0.1)
  expect_equal(results$estimate[2], 0, tolerance = 0.1)
  expect_equal(results$estimate[3], 0, tolerance = 0.1)
})
