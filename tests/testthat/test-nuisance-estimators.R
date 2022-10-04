library(here)
source(here("tests", "testthat", "testing-utils.R"))

test_that(
  paste("fit_prop_score() outputs a vector of estimates that is similar to",
        "ground truth when a validation dataset isn't provided"),
  {

    library(sl3)

    # generate data
    set.seed(5123)
    dt <- generate_test_data(n_obs = 200)

    # fit the propensity score
    fit <- fit_prop_score(train_data = dt,
                          valid_data = NULL,
                          learners = sl3::Lrnr_glm_fast$new(),
                          exposure = "a",
                          confounders = c("w_1", "w_2", "w_3"))

    # make sure that the MSE is approximately zero
    expect_equal(mean((fit$estimates - dt$prop_score)^2), 0, tolerance = 0.01)

})

test_that(
  "fit_prop_score() outputs the sl3 lrnr object fit on the training data",
{

    library(sl3)

    # generate data
    set.seed(5123)
    dt <- generate_test_data(n_obs = 200)

    # fit the propensity score
    fit <- fit_prop_score(train_data = dt,
                          valid_data = NULL,
                          learners = sl3::Lrnr_glm_fast$new(),
                          exposure = "a",
                          confounders = c("w_1", "w_2", "w_3"))

    expect_equal(class(fit$fit), c("Lrnr_glm_fast", "Lrnr_base", "R6"))

})

test_that(
  paste("fit_prop_score() outputs a vector of estimates that is similar to",
        "ground truth when a validation set is provided"),
  {

    library(sl3)

    # generate data
    set.seed(5123)
    train_dt <- generate_test_data(n_obs = 400)
    valid_dt <- generate_test_data(n_obs = 200)

    # fit the propensity score
    fit <- fit_prop_score(train_data = train_dt,
                          valid_data = valid_dt,
                          learners = sl3::Lrnr_glm_fast$new(),
                          exposure = "a",
                          confounders = c("w_1", "w_2", "w_3"))

    # make sure that the MSE is approximately zero
    expect_equal(mean((fit$estimates - valid_dt$prop_score)^2), 0,
                 tolerance = 0.01)

})

test_that(
  paste("fit_cond_est() returns a vector of estimates that is similar to",
        "ground truth when no validation set is provided"),
  {
    library(sl3)

    # generate data
    set.seed(12312)
    dt <- generate_test_data(n_obs = 200)

    # fit the propensity score
    fit <- fit_cond_outcome(train_data = dt,
                          valid_data = NULL,
                          learners = sl3::Lrnr_glm_fast$new(),
                          exposure = "a",
                          confounders = c("w_1", "w_2", "w_3"),
                          outcome = "y")

    # make sure that the MSE is approximately zero
    expect_equal(mean((fit$estimates - dt$y)^2), 0, tolerance = 0.01)
})

test_that(
  paste("fit_cond_est() returns a vector of estimates that is similar to",
        "ground truth when a validation set is provided"),
  {
    library(sl3)

    # generate data
    set.seed(62412)
    train_dt <- generate_test_data(n_obs = 400)
    valid_dt <- generate_test_data(n_obs = 200)

    # fit the propensity score
    fit <- fit_cond_outcome(train_data = train_dt,
                          valid_data = valid_dt,
                          learners = sl3::Lrnr_glm_fast$new(),
                          exposure = "a",
                          confounders = c("w_1", "w_2", "w_3"),
                          outcome = "y")

    # make sure that the MSE is approximately zero
    expect_equal(mean((fit$estimates - valid_dt$y)^2), 0, tolerance = 0.01)
})

test_that(
  "fit_cond_outcome() outputs the sl3 lrnr object fit on the training data",
{

    library(sl3)

    # generate data
    set.seed(5123)
    dt <- generate_test_data(n_obs = 200)

    # fit the propensity score
    fit <- fit_cond_outcome(train_data = dt,
                            valid_data = NULL,
                            learners = sl3::Lrnr_glm_fast$new(),
                            exposure = "a",
                            confounders = c("w_1", "w_2", "w_3"),
                            outcome = "y")

    expect_equal(class(fit$fit), c("Lrnr_glm_fast", "Lrnr_base", "R6"))

})
