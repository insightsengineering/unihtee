test_that("plot_temvip() returns ggplots", {

  # load required packages
  library(MASS)
  library(data.table)
  library(sl3)

  set.seed(843134)

  # create the dataset
  n_obs <- 1000
  w <- mvrnorm(n = n_obs, mu = rep(0, 10), Sigma = diag(10))
  w[, 4] <- ifelse(w[, 4] < 0, 0, 1)
  confounder_names <- paste0("w_", seq_len(10))
  colnames(w) <- confounder_names
  a <- rbinom(n = n_obs, size = 1, prob = plogis(w[, 1] + w[, 2]))
  y <- rnorm(n = n_obs, mean = w[, 1] + w[, 2] + a * w[, 3] - 2 * a * w[, 4])
  dt <- as.data.table(cbind(w, a, y))

  # targeted maximum likelihood estimates and testing procedure
  unihtee_output <- unihtee(
    data = dt,
    confounders = confounder_names,
    modifiers = confounder_names,
    exposure = "a",
    outcome = "y",
    outcome_type = "continuous",
    effect = "absolute",
    estimator = "tmle"
  )

  # plot the TEM-VIP estimate for w_3 and w_4
  expect_no_error(plot_temvip(unihtee_output, modifier_name = "w_3", FALSE))
  expect_no_error(plot_temvip(unihtee_output, modifier_name = "w_4", FALSE))


})
