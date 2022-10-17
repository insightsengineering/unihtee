source("testing-utils.R")

test_that(
  paste(
    "tte_data_melt() melts wide time-to-event data with the correct number",
    "of repeated measures"
  ),
  {

    # generate the data
    n_obs <- 10
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
      propensity_score_values = NULL
    )

    # number of observations that should be contained in the melted data
    unique_times <- unique(dt[time <= time_cutoff, time], time_cutoff)
    reps_per_obs <- sapply(
      seq_len(n_obs),
      function(obs) sum(dt$time[obs] >= unique_times)
    )
    total_reps <- sum(reps_per_obs)

    # assessment
    expect_equal(nrow(melted_dt), total_reps)

  }
)
