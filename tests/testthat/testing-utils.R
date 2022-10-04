# load required libraries
library(data.table)

# randomly generates a basic dataset for testing
generate_test_data <- function(n_obs = 200) {

  # confounders
  w_1 <- rnorm(n = n_obs)
  w_2 <- rnorm(n = n_obs)
  w_3 <- rnorm(n = n_obs)

  # exposure
  prop_score <- plogis(w_1 + w_2)
  a <- rbinom(n = n_obs, size = 1, prob = prop_score)

  # outcome
  y <- rnorm(n = n_obs, mean = w_1 + w_2 + w_3 + a, sd = 0.1)

  # assemble the data.table
  dt <- data.table(
    w_1 = w_1,
    w_2 = w_2,
    w_3 = w_3,
    prop_score = prop_score,
    a = a,
    y = y
  )

  return(dt)
}
