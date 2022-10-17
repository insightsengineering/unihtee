utils::globalVariables(c(
  "modifier", "estimate", "var", "se",
  "z", "p_value", "ci_lower", "ci_upper", "p_value_fdr"
))

#' @title Test Hypotheses
#'
#' @description \code{test_hypotheses()} is a convenience function for testing
#'   null hypotheses.
#'
#' @param n_obs A \code{numeric} representing the number of observations in the
#'   dataset.
#' @param estimates A one-row \code{data.table} of estimates for each potential
#'   treatment effect modifier.
#' @param var_estimates A one-row \code{data.table} of estimator variances for
#'   each potential treatment effect modifier.
#' @param rescale_factor A \code{numeric} factor by which to re-scale the
#'   estimates and efficient influence functions. Used when the outcome
#'   variable is continuous but not bounded between 0 and 1.
#'
#' @return A \code{data.table} containing the hypothesis testing results.
#'
#' @importFrom data.table melt ":="
#' @importFrom stats pnorm p.adjust
#'
#' @keywords internal
test_hypotheses <- function(n_obs, estimates, var_estimates, rescale_factor) {

  # pivot estimates dt to long format
  # NOTE: How do I get rid of these annoying warnings?
  estimates <- estimates * rescale_factor
  estimates_l <- suppressWarnings(data.table::melt(
    data = estimates,
    variable.name = "modifier",
    value.name = "estimate"
  ))

  # pivot variance dt to long format
  var_estimates <- var_estimates * rescale_factor^2
  var_estimates_l <- suppressWarnings(data.table::melt(
    data = var_estimates,
    variable.name = "modifier",
    value.name = "var"
  ))

  # mege dts
  test_dt <- estimates_l[var_estimates_l, on = "modifier"]

  # compute the standard errors, Z-score, p values and adjusted p values
  test_dt[, se := sqrt(var / n_obs), by = modifier]
  test_dt[, z := estimate / se, by = modifier]
  test_dt[, p_value := (2 * min(stats::pnorm(z), 1 - stats::pnorm(z))),
    by = modifier
  ]
  test_dt[, ci_lower := estimate - 1.96 * se, by = modifier]
  test_dt[, ci_upper := estimate + 1.96 * se, by = modifier]
  test_dt[, p_value_fdr := stats::p.adjust(p_value, method = "BH")]
  test_dt[, var := NULL]

  return(test_dt)
}
