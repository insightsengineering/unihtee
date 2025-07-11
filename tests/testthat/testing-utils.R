# load required libraries
library(data.table)

# randomly generates a basic dataset for testing
generate_test_data <- function(
  n_obs = 200,
  outcome_type = "continuous",
  centered_modifiers = TRUE
) {

  # confounders
  if (centered_modifiers) {
    w_1 <- rnorm(n = n_obs)
    w_2 <- rnorm(n = n_obs)
    w_3 <- rnorm(n = n_obs)
  } else {
    w_1 <- rnorm(n = n_obs, mean = 1)
    w_2 <- rnorm(n = n_obs, mean = -1)
    w_3 <- rbinom(n = n_obs, size = 1, prob = 0.5)
  }

  # exposure
  prop_score <- plogis(w_1 + w_2)
  a <- rbinom(n_obs, 1, prop_score)

  # DGP type
  if (outcome_type %in% c("continuous", "binary")) {
    if (outcome_type == "continuous") {
      # (potential) outcomes
      y_1 <- rnorm(n = n_obs, mean = w_1 + w_2 + w_3 + 1, sd = 0.1)
      y_0 <- rnorm(n = n_obs, mean = w_1 + w_2, sd = 0.1)
    } else if (outcome_type == "binary") {
      resp_prob_1 <- plogis(-1 + 2 * w_3)
      resp_prob_0 <- plogis(1 - 2 * w_3)
      y_1 <- rbinom(n_obs, 1, resp_prob_1)
      y_0 <- rbinom(n_obs, 1, resp_prob_0)
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

  } else if (outcome_type %in% c("time-to-event", "time-to-event RR")) {

    # define hazard functions
    cond_cens_hazard <- function(time, exposure, w_1, w_2, w_3) 0.05
    if (outcome_type == "time-to-event RR") {
      cond_surv_hazard <- function(time, exposure, w_1, w_2, w_3) {
        (time < 9) / (2 + exp(2 * exposure * w_1)) + (time == 9)
      }
    } else {
      # very prone to survival probabilities near 0 or 1
      cond_surv_hazard <- function(time, exposure, w_1, w_2, w_3) {
        (time < 9) / (1 + exp(2 + 3 * exposure * w_1)) + (time == 9)
      }
    }

    # generate the failure events for t = 1 to 9
    failure_time_sim <- function(exposure) {
      sapply(
        seq_len(n_obs),
        function(obs) {
          failure_time <- NA
          for (t in 1:9) {
            prob <- cond_surv_hazard(t, exposure, w_1[obs], w_2[obs], w_3[obs])
            status <- rbinom(1, 1, prob)
            if (status == 1) {
              failure_time <- t
              break
            }
          }
          return(failure_time)
        }
      )
    }
    failure_time_1 <- failure_time_sim(1)
    failure_time_0 <- failure_time_sim(0)


    # generate the censoring events for t = 1 to 9
    censor_time_sim <- function(exposure) {
      sapply(
        seq_len(n_obs),
        function(obs) {
          censor_time <- NA
          for (t in 1:9) {
            prob <- cond_cens_hazard(t, exposure, w_1[obs], w_2[obs], w_3[obs])
            status <- rbinom(1, 1, prob)
            if (status == 1) {
              censor_time <- t
              break
            }
          }
          if (is.na(censor_time)) censor_time <- 10
          return(censor_time)
        }
      )
    }
    censor_time_1 <- censor_time_sim(1)
    censor_time_0 <- censor_time_sim(0)

    # compile the failure and censoring times
    failure_time <- sapply(
      seq_len(n_obs),
      function(obs) {
        if (a[obs] == 1) failure_time_1[obs] else failure_time_0[obs]
      }
    )
    censor_time <- sapply(
      seq_len(n_obs),
      function(obs) {
        if (a[obs] == 1) censor_time_1[obs] else censor_time_0[obs]
      }
    )

    # assess the observed time-to-event and censoring indicator
    time <- sapply(
      seq_len(n_obs),
      function(obs) {
        if (censor_time[obs] < failure_time[obs]) {
          censor_time[obs]
        } else {
          failure_time[obs]
        }
      }
    )
    censoring <- sapply(
      seq_len(n_obs),
      function(obs) if (time[obs] == censor_time[obs]) 1 else 0
    )

    # assemble the data.table
    dt <- data.table(
      w_1 = w_1,
      w_2 = w_2,
      w_3 = w_3,
      prop_score = prop_score,
      a = a,
      time = time,
      censoring = censoring
    )


  }

  return(dt)
}
