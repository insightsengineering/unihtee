#' @title Time-to-Event Data Melt
#'
#' @description \code{tte_data_melt()} turns wide-format \code{data.table}
#'   objects into long-format \code{data.table} objects.
#'
#' @param data A \code{data.table} containing the observed data.
#' @param confounders A \code{character} vector of column names corresponding to
#'   baseline covariates.
#' @param exposure A \code{character} corresponding to the exposure variable.
#' @param outcome A \code{character} corresponding to the outcome variable.
#' @param censoring A \code{character} indicating the right censogin indicator
#'   variable.
#' @param time_cutoff A \code{numeric} representing the timepoint at which to
#'   evaluate the time-to-event parameter.
#' @param prop_score_values A \code{character} corresponding to the (known)
#'   propensity score values for each observation in \code{data}.
#'
#' @return A long, time-to-event version of the original \code{data.table}
#'   object. That is, observations are repeated at timepoints between the
#'   earliest time in the data and the \code{time_cutoff} argument.
#'
#' @importFrom data.table .N rbindlist
#'
#' @keywords internal
#'
tte_data_melt <- function(data,
                          confounders,
                          exposure,
                          outcome,
                          censoring,
                          time_cutoff,
                          prop_score_values) {

  ## remove unecessary features from the data
  to_keep <- c(confounders, exposure, outcome, censoring, prop_score_values)
  data <- data[, ..to_keep]

  ## find the number of unique times up to and including the time_cutoff
  unique_times <- sort(
    unique(c(data[get(outcome) < time_cutoff, get(outcome)], time_cutoff))
  )

  ## expand each observations data individually
  obs_ls <- lapply(
    seq_len(nrow(data)),
    function(obs_idx) {

      ## extract the observation
      obs_data <- data[obs_idx, ]

      ## check if it's censored
      censored <- obs_data[[censoring]]
      obs_data[[censoring]] <- 0

      ## set the failure indicator to 0
      obs_data$failure <- 0

      ## store event time for later
      event_time <- obs_data[[outcome]]

      ## count the number of repititions
      ## NOTE: add a rep at the time_cutoff if an event occurs after it
      num_reps <- sum(obs_data[[outcome]] >= unique_times)

      ## repeat the observation num_reps times
      obs_data <- obs_data[rep(1:.N, num_reps)]

      ## add the pseudo observation times
      obs_data[[outcome]] <- unique_times[1:num_reps]

      ## indicate if censored at the last measured timepoint
      if (censored && event_time <= time_cutoff) {
        obs_data[num_reps, censoring] <- 1
      }

      ## indicate failure at the last measured timepoint
      if (censored == 0 && event_time <= time_cutoff) {
        failure <- "failure"
        obs_data[num_reps, failure] <- 1
      }

      ## should the repeated measure be kept when estimating nuisance
      ## parameters?
      obs_data$keep <- 1

      ## if observation does not make it to designated timepoint, add dummy
      ## observations used for computing at that timepoint anyways
      ## NOTE: these repeated measures are dropped during nuisance parameter
      ## estimation, but used for EIF calculation
      if (num_reps < length(unique_times)) {
        num_add_reps <- length(unique_times) - num_reps
        add_obs_data <- data[obs_idx]
        add_obs_data$failure <- 0
        add_obs_data$keep <- 0
        add_obs_data <- add_obs_data[rep(1:.N, num_add_reps)]
        obs_data <- data.table::rbindlist(list(obs_data, add_obs_data))
        obs_data[[outcome]] <- unique_times
      }

      return(obs_data)
    }
  )

  ## combine into a single long data.table
  data.table::rbindlist(obs_ls, idcol = "id")

}
