source("testing-utils.R")

test_that(
  paste(
    "tte_data_melt() melts wide time-to-event data with the correct number",
    "of repeated measures"
  ),
  {

    set.seed(72345)

    # generate the data
    n_obs <- 15
    dt <- generate_test_data(n_obs = n_obs, outcome_type = "time-to-event")

    # melt the data appropriately
    time_cutoff <- 5
    melted_dt <- tte_data_melt(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = time_cutoff,
      prop_score_values = NULL
    )

    # number of observations that should be contained in the melted data
    unique_times <- unique(c(dt[time <= time_cutoff, time], time_cutoff))
    reps_per_obs <- sapply(
      seq_len(n_obs),
      function(obs) sum(dt$time[obs] >= unique_times)
    )
    total_reps <- sum(reps_per_obs)

    expect_equal(sum(melted_dt$keep), total_reps)
    expect_equal(nrow(melted_dt), length(unique_times) * n_obs)

  }
)

test_that(
  paste(
    "tte_data_melt() outputs a data.table with columns for confounders,",
    "the exposure, the outcome, the censoring indicator, and if necessary",
    "the known propensity scores"
  ),
  {

    set.seed(15412)

    # generate the data
    n_obs <- 15
    dt <- generate_test_data(n_obs = n_obs, outcome_type = "time-to-event")
    dt$ps <- rep(0.5, 15)

    # melt the data appropriately
    time_cutoff <- 5
    melted_dt <- tte_data_melt(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = time_cutoff,
      prop_score_values = "ps"
    )

    expect_equal(
      colnames(melted_dt),
      c("id", "w_1", "w_2", "w_3", "a", "time", "censoring", "ps", "failure",
        "keep")
    )

    melted_dt <- tte_data_melt(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = time_cutoff,
      prop_score_values = NULL
    )

    expect_equal(
      colnames(melted_dt),
      c("id", "w_1", "w_2", "w_3", "a", "time", "censoring", "failure", "keep")
    )
  }
)
