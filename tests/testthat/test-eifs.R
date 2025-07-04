source("testing-utils.R")

test_that(
  paste(
    "uncentered_eif() minus true absolute effect parameter value has a mean",
    "of zero when propensity scores aren't known"
  ),
  {
    library(sl3)

    # generate data
    set.seed(5123)
    dt <- generate_test_data(n_obs = 100000)

    # fit the propensity score
    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(
        formula = "~ w_1 + w_2 + w_3"
      ),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected cond outcome
    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(
        formula = "~ w_1 * a + w_2 * a + w_3 * a"
      ),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    # compute the uncentered eif
    eif <- uncentered_eif(
      data = dt,
      effect = "absolute",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = prop_score_fit,
      prop_score_values = NULL,
      cond_outcome_fit = cond_outcome_fit,
      failure_hazard_fit = NULL,
      censoring_hazard_fit = NULL
    )

    # note that the true parameter values are equal to 0,1 for W_1, W_3
    expect_equal(eif[, sapply(.SD, mean)], c("w_1" = 0, "w_3" = 1),
      tolerance = 0.005
    )
  }
)

test_that(
  paste(
    "uncentered_eif() minus true absolute effect parameter value has a mean",
    "of zero when propensity scores are known"
  ),
  {
    library(sl3)

    # generate data
    set.seed(12723)
    dt <- generate_test_data(n_obs = 5000)

    # fit the expected cond outcome
    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glmnet$new(
        formula = "~ a * w_1 + a * w_2 + a * w_3"
      ),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    # compute the uncentered eif
    eif <- uncentered_eif(
      data = dt,
      effect = "absolute",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = NULL,
      cond_outcome_fit = cond_outcome_fit,
      prop_score_values = "prop_score",
      failure_hazard_fit = NULL,
      censoring_hazard_fit = NULL
    )

    # note that the true parameter values are equal to 0,1 for W_1, W_3
    expect_equal(eif[, sapply(.SD, mean)], c("w_1" = 0, "w_3" = 1),
      tolerance = 0.1
    )
  }
)

test_that(
  paste(
    "uncentered_eif() minus true relative effect parameter value has a mean",
    "of zero when propensity scores aren't known"
  ),
  {
    library(sl3)

    # generate data
    set.seed(124151)
    dt <- generate_test_data(n_obs = 20000, outcome_type = "binary")

    # fit the propensity score
    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected cond outcome
    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_xgboost$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    # compute the uncentered eif
    eif <- uncentered_eif(
      data = dt,
      effect = "relative",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = prop_score_fit,
      prop_score_values = NULL,
      cond_outcome_fit = cond_outcome_fit,
      failure_hazard_fit = NULL,
      censoring_hazard_fit = NULL
    )

    # note that the true parameter values are approx equal to 0, 2.0 for W_1,
    # W_3
    expect_equal(eif[, sapply(.SD, mean)], c("w_1" = 0, "w_3" = 2.0),
      tolerance = 0.1
    )
  }
)

test_that(
  paste(
    "uncentered_eif() minus true relative effect parameter value has a mean",
    "of zero when propensity scores are known"
  ),
  {
    library(sl3)

    # generate data
    set.seed(82342)
    dt <- generate_test_data(n_obs = 20000, outcome_type = "binary")

    # fit the expected cond outcome
    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_xgboost$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    # compute the uncentered eif
    eif <- uncentered_eif(
      data = dt,
      effect = "relative",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = NULL,
      cond_outcome_fit = cond_outcome_fit,
      prop_score_values = "prop_score",
      failure_hazard_fit = NULL,
      censoring_hazard_fit = NULL
    )

    # note that the true parameter values are equal approx to 0, 2.0 for W_1,
    # W_3
    expect_equal(eif[, sapply(.SD, mean)], c("w_1" = 0, "w_3" = 2.0),
      tolerance = 0.1
    )
  }
)

test_that(
  paste(
    "uncentered_eif() minus true absolute effect parameter value has a mean",
    "of zero for the RD TTE TEM VIP"
  ),
  {
    library(sl3)

    # generate data
    set.seed(16234)
    dt <- generate_test_data(n_obs = 10000, outcome_type = "time-to-event")
    long_dt <- tte_data_melt(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = 5,
      prop_score_values = NULL
    )

    # fit the propensity score
    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = long_dt,
      learners = sl3::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected failure hazard
    fail_fit <- fit_failure_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3:::Lrnr_xgboost$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      times = "time"
    )

    # fit the expected failure hazard
    cens_fit <- fit_censoring_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3:::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      times = "time",
      censoring = "censoring"
    )

    # compute the uncentered eif
    eif <- uncentered_eif(
      data = long_dt,
      effect = "absolute",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = prop_score_fit,
      cond_outcome_fit = NULL,
      prop_score_values = NULL,
      failure_hazard_fit = fail_fit,
      censoring_hazard_fit = cens_fit
    )

    # get approximations to the true parameter values

    ## # generate the data
    ## set.seed(4514)
    ## dt <- generate_test_data(n_obs = 100000, outcome_type = "time-to-event")

    # make it long
    long_dt <- lapply(
      seq_len(nrow(dt)),
      function(obs) {
        obs_dt <- dt[obs]
        obs_dt <- obs_dt[rep(1:.N, 5)]
        obs_dt$time <- seq_len(5)
        obs_dt
      }
    )
    long_dt <- rbindlist(long_dt, idcol = "id")

    # compute the true failure hazard at each timepoint under each condition
    cond_surv_hazard <- function(time, exposure, w_1, w_2, w_3) {
      (time < 9) / (1 + exp(2 + 3 * exposure * w_1)) + (time == 9)
    }
    exp_truth <-
        cond_surv_hazard(long_dt$time, 1, long_dt$w_1, long_dt$w_2, long_dt$w_3)
    noexp_truth <-
        cond_surv_hazard(long_dt$time, 0, long_dt$w_1, long_dt$w_2, long_dt$w_3)

    # compute the survival probability at each time
    long_dt$true_haz_exp <- exp_truth
    long_dt$true_haz_noexp <- noexp_truth
    long_dt[, surv_exp := cumprod(1 - true_haz_exp), by = "id"]
    long_dt[, surv_noexp := cumprod(1 - true_haz_noexp), by = "id"]

    # compute the restricted mean survival time differences
    long_dt[, rmst := cumsum(surv_exp - surv_noexp), by = "id"]

    # retain only the rmst differences at time_cutoff, and compute parameters
    res_dt <- long_dt[time == 5]
    w_1_param <- cov(res_dt$w_1, res_dt$rmst) / var(res_dt$w_1) # 1.721
    w_3_param <- cov(res_dt$w_3, res_dt$rmst) / var(res_dt$w_3) # -0.005

    expect_equal(mean(eif$w_1), w_1_param, tolerance = 0.1)
    expect_equal(mean(eif$w_3), w_3_param, tolerance = 0.1)
  }
)

test_that(
  paste(
    "uncentered_eif() minus true relative effect parameter value has a mean",
    "of zero for the RR TTE TEM VIP"
  ),
  {
    library(sl3)

    # generate data
    set.seed(42)
    dt <- generate_test_data(n_obs = 10000, outcome_type = "time-to-event RR")
    long_dt <- tte_data_melt(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      censoring = "censoring",
      time_cutoff = 3,
      prop_score_values = "prop_score"
    )

    # fit the expected failure hazard
    fail_fit <- fit_failure_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3:::Lrnr_xgboost$new(),
      exposure = "a",
      times = "time",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected failure hazard
    cens_fit <- fit_censoring_hazard(
      train_data = long_dt,
      valid_data = NULL,
      learners = sl3:::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      times = "time",
      censoring = "censoring"
    )

    eif <- uncentered_eif(
      data = long_dt,
      effect = "relative",
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "time",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = NULL,
      cond_outcome_fit = NULL,
      prop_score_values = "prop_score",
      failure_hazard_fit = fail_fit,
      censoring_hazard_fit = cens_fit
    )

    ## # get approximations to the true parameter values

    ## #  manually generate long data
    ## long_dt <- lapply(
    ##   seq_len(nrow(dt)),
    ##   function(obs) {
    ##     obs_dt <- dt[obs]
    ##     obs_dt <- obs_dt[rep(1:.N, 3)]
    ##     obs_dt$time <- seq_len(3)
    ##     obs_dt
    ##   }
    ## )
    ## long_dt <- rbindlist(long_dt, idcol = "id")

    ## # compute the true failure hazard at each timepoint under each condition
    ## cond_surv_hazard <- function(time, exposure, w_1, w_2, w_3) {
    ##   (time <= 9) / (2 + exp(2 * exposure * w_1))
    ## }
    ## exp_truth <-
    ##     cond_surv_hazard(long_dt$time, 1, long_dt$w_1, long_dt$w_2, long_dt$w_3)
    ## noexp_truth <-
    ##     cond_surv_hazard(long_dt$time, 0, long_dt$w_1, long_dt$w_2, long_dt$w_3)

    ## # compute the survival probability at each time under each condition
    ## long_dt$true_haz_exp <- exp_truth
    ## long_dt$true_haz_noexp <- noexp_truth
    ## long_dt[, surv_exp := cumprod(1 - true_haz_exp), by = "id"]
    ## long_dt[, surv_noexp := cumprod(1 - true_haz_noexp), by = "id"]

    ## # retain only the log ratio at time_cutoff, and compute parameters
    ## res_dt <- long_dt[time == 3]
    ## res_dt[, log_ratio := log(surv_exp / surv_noexp), by = "id"]
    ## w_1_param <- cov(res_dt$w_1, res_dt$log_ratio) / var(res_dt$w_1) # 0.6
    ## w_3_param <- cov(res_dt$w_3, res_dt$log_ratio) / var(res_dt$w_3) # 0.0

    expect_equal(mean(eif$w_1) - 0.62, 0, tolerance = 0.1)
    expect_equal(mean(eif$w_3), 0, tolerance = 0.1)
  }
)

test_that(
  paste(
    "uncentered_eif() produces standard errors near zero in large datasets"
  ),
  {
    library(sl3)

    # generate data
    set.seed(5123)
    sample_size <- 100000
    dt <- generate_test_data(n_obs = sample_size, centered_modifiers = TRUE)

    # fit the propensity score
    prop_score_fit <- fit_prop_score(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::Lrnr_glm_fast$new(),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3")
    )

    # fit the expected cond outcome
    cond_outcome_fit <- fit_cond_outcome(
      train_data = dt,
      valid_data = NULL,
      learners = sl3::make_learner(
        sl3::Pipeline,
        sl3::Lrnr_define_interactions$new(
          list(c("w_1", "a"), c("w_2", "a"), c("w_3", "a"))
        ),
        sl3::Lrnr_glm_fast$new()
      ),
      exposure = "a",
      confounders = c("w_1", "w_2", "w_3"),
      outcome = "y"
    )

    # estimate the ATE
    ate_estimate <- one_step_ate_estimator(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      prop_score_fit = prop_score_fit,
      prop_score_values = NULL,
      cond_outcome_fit = cond_outcome_fit
    )

    # compute the eif
    eif <- centered_eif(
      data = dt,
      confounders = c("w_1", "w_2", "w_3"),
      exposure = "a",
      outcome = "y",
      modifiers = c("w_1", "w_3"),
      prop_score_fit = prop_score_fit,
      prop_score_values = NULL,
      cond_outcome_fit = cond_outcome_fit,
      ace_estimate = ate_estimate
    )

    expect_equal(sqrt(eif[, sapply(.SD, var)] / sample_size),
                 c("w_1" = 0, "w_3" = 0),
                 tolerance = 0.01
    )
  }
)
