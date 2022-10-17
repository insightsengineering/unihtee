# load required libraries
library(data.table)

# randomly generates a basic dataset for testing
generate_test_data <- function(n_obs = 200, outcome_type = "continuous") {

  # confounders
  w_1 <- rnorm(n = n_obs)
  w_2 <- rnorm(n = n_obs)
  w_3 <- rnorm(n = n_obs)

  # exposure
  prop_score <- plogis(w_1 + w_2)
  a <- sapply(prop_score, function(x) rbinom(1, 1, prob = x))

  if (outcome_type == "continuous") {
    # (potential) outcomes
    y_1 <- rnorm(n = n_obs, mean = w_1 + w_2 + w_3 + 1, sd = 0.1)
    y_0 <- rnorm(n = n_obs, mean = w_1 + w_2, sd = 0.1)
  } else if (outcome_type == "binary") {
    resp_prob_1 <- plogis(2 + 10 * w_3)
    resp_prob_0 <- plogis(1)
    y_1 <- sapply(resp_prob_1, function(x) rbinom(1, 1, x))
    y_0 <- sapply(resp_prob_0, function(x) rbinom(1, 1, x))
  }

  # outcome
  y <- a * y_1 + (1 - a) * y_0

  # assemble the data.table
  dt <- data.table(
    w_1 = w_1,
    w_2 = w_2,
    w_3 = w_3,
    prop_score = prop_score,
    a = a,
    y_1 = y_1,
    y_0 = y_0,
    y = y
  )

  return(dt)
}
