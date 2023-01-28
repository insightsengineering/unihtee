source("testing-utils.R")

test_that(
  paste(
    "unihtee() uncovers treatment effect modifiers on the absolute effect",
    "when outcomes are continuous with the one-step estimator"
  ),
  {

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
      effect = "absolute",
      estimator = "onestep"
    )

    # ensure that the adjusted p-value of w_3 is less than 0.05, and those of
    # w_1 and w_2 are above 0.05
    expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))

    # check that the estimates are reasonably close to the ground truth
    expect_equal(results$estimate[1], 1, tolerance = 0.1)
    expect_equal(results$estimate[2], 0, tolerance = 0.1)
    expect_equal(results$estimate[3], 0, tolerance = 0.1)
  }
)


test_that(
  paste(
    "unihtee() uncovers treatment effect modifiers on the absolute effect",
    "when outcomes are continuous with the cross-fitted one-step estimator"
  ),
  {

    # generate data
    set.seed(8234)
    dt <- generate_test_data(n_obs = 500)

    # apply unihtee
    results <- unihtee(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      modifiers = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      outcome_type = "continuous",
      effect = "absolute",
      estimator = "onestep",
      cross_fit = TRUE,
      cond_outcome_estimator = sl3::Lrnr_xgboost$new()
    )

    # ensure that the adjusted p-value of w_3 is less than 0.05, and those of
    # w_1 and w_2 are above 0.05
    expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))
  }
)


test_that(
  paste(
    "unihtee() uncovers treatment effect modifiers on the absolute effect",
    "when outcomes are continuous with the TML estimator"
  ),
  {

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
      effect = "absolute",
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
  }
)

test_that(
  paste(
    "unihtee() uncovers treatment effect modifiers on the absolute effect",
    "when outcomes are continuous with the cross-fitted TML estimator"
  ),
  {

    # generate data
    set.seed(12345)
    dt <- generate_test_data(n_obs = 200)

    # apply unihtee
    results <- unihtee(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      modifiers = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      outcome_type = "continuous",
      effect = "absolute",
      estimator = "tmle",
      cross_fit = TRUE,
      cond_outcome_estimator = sl3::Lrnr_glmnet$new(
        formula = "~ a * w_1 + a * w_2 + a * w_3"
      )
    )

    # ensure that the adjusted p-value of w_3 is less than 0.05, and those of
    # w_1 and w_2 are above 0.05
    expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))
  }
)

test_that(
  paste(
    "unihtee() uncovers treatment effect modifiers on the absolute effect",
    "when outcomes are binary with the one-step estimator"
  ),
  {

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
      effect = "absolute",
      estimator = "onestep"
    )

    # ensure that the adjusted p-value of w_3 is less than 0.05, and those of
    # w_1 and w_2 are above 0.05
    expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))

    # check that the estimates are reasonably close to the ground truth
    expect_equal(results$estimate[1], 0.57, tolerance = 0.1)
    expect_equal(results$estimate[2], 0, tolerance = 0.1)
    expect_equal(results$estimate[3], 0, tolerance = 0.1)
  }
)

test_that(
  paste(
    "unihtee() uncovers treatment effect modifiers on the absolute effect",
    "when outcomes are binary with the TML estimator"
  ),
  {

    # generate data
    set.seed(823421)
    dt <- generate_test_data(n_obs = 5000, outcome_type = "binary")

    # apply unihtee
    results <- unihtee(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      modifiers = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      outcome_type = "binary",
      effect = "absolute",
      estimator = "tmle",
      cond_outcome_estimator = sl3::Lrnr_xgboost$new(),
      prop_score_estimator = sl3::Lrnr_xgboost$new()
    )

    # ensure that the adjusted p-value of w_1 is less than 0.05, and those of
    # w_2 and w_2 are above 0.05
    expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))

    # check that the estimates are reasonably close to the ground truth
    expect_equal(results$estimate[1], 0.57, tolerance = 0.1)
    expect_equal(results$estimate[2], 0, tolerance = 0.1)
    expect_equal(results$estimate[3], 0, tolerance = 0.1)
  }
)

test_that(
  paste(
    "unihtee() uncovers treatment effect modifiers on the relative effect",
    "when outcomes are binary with the one-step estimator"
  ),
  {

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
      effect = "relative",
      estimator = "onestep",
      cond_outcome_estimator = sl3:::Lrnr_xgboost$new()
    )

    # ensure that the adjusted p-value of w_3 is less than 0.05, and those of
    # w_1 and w_2 are above 0.05
    expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))

    # check that the estimates are reasonably close to the ground truth
    expect_equal(results$estimate[1], 2.0, tolerance = 0.1)
    expect_equal(results$estimate[2], 0, tolerance = 0.1)
    expect_equal(results$estimate[3], 0, tolerance = 0.1)
  }
)

test_that(
  paste(
    "unihtee() uncovers treatment effect modifiers on the relative effect",
    "when outcomes are binary with the cross-fitted one-step estimator"
  ),
  {

    # generate data
    set.seed(90234)
    dt <- generate_test_data(n_obs = 500, outcome_type = "binary")

    # apply unihtee
    results <- unihtee(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      modifiers = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      outcome_type = "binary",
      effect = "relative",
      estimator = "onestep",
      cross_fit = TRUE
    )

    # ensure that the adjusted p-value of w_3 is less than 0.05, and those of
    # w_1 and w_2 are above 0.05
    expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))
  }
)

test_that(
  paste(
    "unihtee() uncovers treatment effect modifiers on the relative effect",
    "when outcomes are binary with the TML estimator"
  ),
  {

    # generate data
    set.seed(5140)
    dt <- generate_test_data(n_obs = 5000, outcome_type = "binary")

    # apply unihtee
    results <- unihtee(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      modifiers = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      outcome_type = "binary",
      effect = "relative",
      estimator = "tmle"
    )

    # ensure that the adjusted p-value of w_3 is less than 0.05, and those of
    # w_1 and w_2 are above 0.05
    expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))

    # check that the estimates are reasonably close to the ground truth
    expect_equal(results$estimate[1], 2.0, tolerance = 0.1)
    expect_equal(results$estimate[2], 0, tolerance = 0.1)
    expect_equal(results$estimate[3], 0, tolerance = 0.1)
  }
)

test_that(
  paste(
    "unihtee() uncovers treatment effect modifiers on the relative effect",
    "when outcomes are binary with the cross-fitted TML estimator"
  ),
  {

    # generate data
    set.seed(5140)
    dt <- generate_test_data(n_obs = 500, outcome_type = "binary")

    # apply unihtee
    results <- unihtee(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      modifiers = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      outcome_type = "binary",
      effect = "relative",
      estimator = "tmle",
      cross_fit = TRUE
    )

    # ensure that the adjusted p-value of w_3 is less than 0.05, and those of
    # w_1 and w_2 are above 0.05
    expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))
  }
)

test_that(
  paste(
    "unihtee() uncovers treatment effect modifiers on the absolute effect",
    "with time-to-event outcomes with the one-step estimator"
  ),
  {

    # generate data
    set.seed(72342)
    dt <- generate_test_data(n_obs = 5000, outcome_type = "time-to-event")

    # apply unihtee
    results <- unihtee(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      modifiers = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = 5,
      outcome_type = "time-to-event",
      effect = "absolute",
      estimator = "onestep",
      prop_score_estimator = sl3::Lrnr_xgboost$new(),
      cond_outcome_estimator = NULL,
      failure_hazard_estimator = sl3::Lrnr_xgboost$new(),
      censoring_hazard_estimator = sl3::Lrnr_xgboost$new()
    )

    # ensure that the adjusted p-value of w_1 is less than 0.05, and those of
    # w_2 and w_2 are above 0.05
    expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))

    # check that the estimates are reasonably close to the ground truth
    expect_equal(results$estimate[1], 1.72, tolerance = 0.1)
    expect_equal(results$estimate[2], 0, tolerance = 0.1)
    expect_equal(results$estimate[3], 0, tolerance = 0.1)
  }
)

test_that(
  paste(
    "unihtee() uncovers treatment effect modifiers on the absolute effect",
    "with time-to-event outcomes with the cross-fitted one-step estimator"
  ),
  {

    # generate data
    set.seed(517)
    dt <- generate_test_data(n_obs = 500, outcome_type = "time-to-event")

    # apply unihtee
    results <- unihtee(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      modifiers = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = 5,
      outcome_type = "time-to-event",
      effect = "absolute",
      estimator = "onestep",
      cross_fit = TRUE,
      prop_score_estimator = sl3::Lrnr_xgboost$new(),
      cond_outcome_estimator = NULL,
      failure_hazard_estimator = sl3::Lrnr_xgboost$new(),
      censoring_hazard_estimator = sl3::Lrnr_xgboost$new()
    )

    # ensure that the adjusted p-value of w_1 is less than 0.05, and those of
    # w_2 and w_2 are above 0.05
    expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))
  }
)

test_that(
  paste(
    "unihtee() uncovers treatment effect modifiers on the absolute effect",
    "with time-to-event outcomes with the TML estimator"
  ),
  {

    # generate data
    set.seed(72341)
    dt <- generate_test_data(n_obs = 5000, outcome_type = "time-to-event")

    # apply unihtee
    results <- unihtee(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      modifiers = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = 5,
      outcome_type = "time-to-event",
      effect = "absolute",
      estimator = "tmle",
      prop_score_estimator = sl3::Lrnr_xgboost$new(),
      cond_outcome_estimator = NULL,
      failure_hazard_estimator = sl3::Lrnr_xgboost$new(),
      censoring_hazard_estimator = sl3::Lrnr_xgboost$new()
    )

    # ensure that the adjusted p-value of w_1 is less than 0.05, and those of
    # w_2 and w_2 are above 0.05
    expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))

    # check that the estimates are reasonably close to the ground truth
    expect_equal(results$estimate[1], 1.72, tolerance = 0.1)
    expect_equal(results$estimate[2], 0, tolerance = 0.1)
    expect_equal(results$estimate[3], 0, tolerance = 0.1)
  }
)

test_that(
  paste(
    "unihtee() uncovers treatment effect modifiers on the absolute effect",
    "with time-to-event outcomes with the cross-fitted TML estimator"
  ),
  {

    # generate data
    set.seed(9951)
    dt <- generate_test_data(n_obs = 500, outcome_type = "time-to-event")

    # apply unihtee
    results <- unihtee(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      modifiers = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = 5,
      outcome_type = "time-to-event",
      effect = "absolute",
      estimator = "tmle",
      cross_fit = TRUE,
      prop_score_estimator = sl3::Lrnr_xgboost$new(),
      cond_outcome_estimator = NULL,
      failure_hazard_estimator = sl3::Lrnr_xgboost$new(),
      censoring_hazard_estimator = sl3::Lrnr_xgboost$new()
    )

    # ensure that the adjusted p-value of w_1 is less than 0.05, and those of
    # w_2 and w_2 are above 0.05
    expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))
  }
)

test_that(
  paste(
    "unihtee() uncovers treatment effect modifiers on the relative effect",
    "with time-to-event outcomes with the cross-fitted TML estimator"
  ),
  {

    library(earth)

    # generate data
    set.seed(6234)
    dt <- generate_test_data(n_obs = 1000, outcome_type = "time-to-event")

    # apply unihtee
    results <- unihtee(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      modifiers = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = 5,
      outcome_type = "time-to-event",
      effect = "relative",
      estimator = "tmle",
      cross_fit = TRUE,
      prop_score_estimator = sl3::Lrnr_glm$new(),
      cond_outcome_estimator = NULL,
      failure_hazard_estimator = sl3::Lrnr_earth$new(
        formula = "~ a * w_1 + w_2 + w_3",
        glm = list(family = "binomial")
      ),
      censoring_hazard_estimator = sl3::Lrnr_glm$new()
    )

    # ensure that the adjusted p-value of w_1 is less than 0.05, and those of
    # w_2 and w_2 are above 0.05
    expect_equal(results$p_value_fdr < 0.05, c(TRUE, FALSE, FALSE))
  }
)
