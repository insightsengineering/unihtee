#' Plot Treatment Effect Modifier Variable Importance Parameter Estimate
#'
#' @description \code{plot.unihtee()} produces a plot depicting the estimated
#' simple linear regression line for the specified \code{modifier} parameter
#' estimated by \code{\link{unihtee}()}.
#'
#' @param unihtee_output A \code{unihtee} class object output by
#'   \code{\link{unihtee}()}.
#' @param data A \code{data.table} containing the data used to produce the
#'   \code{unihtee_output} object.
#' @param modifier A \code{character} specifying the treatment effect modifier
#' variable importance parameter estimate to plot.
#'
#' @importFrom ggplot2 ggplot aes
#'
#' @returns A \code{\link[ggplot2]{ggplot2}} object.
#' @export
plot <- function(unihtee_output, data, modifier_name) {
  UseMethod("plot")
}

#' @export
plot.unihtee <- function(unihtee_output, data, modifier_name) {

  # extract the ACE estimate
  ace_estimate <- unihtee_output$ace_estimate

  # extract the TEM-VIP inference
  temvip_inference <- unihtee_output$temvip_inference_tbl[
    modifier == modifier_name, ]
  temvip_est <- temvip_inference$estimate
  temvip_ci_lower <- temvip_inference$ci_lower
  temvip_ci_upper <- temvip_inference$ci_upper

  # extract the modifier values from the data
  modifier_vals <- data[[modifier_name]]

  # compute the mean modifier value
  modifier_mean <- mean(modifier_vals)

  # create the data.frame for plotting
  plotting_tbl <- data.frame(
    x = modifier_vals,
    y = ace_estimate + temvip_est * (modifier_vals - modifier_mean)
  )

  # prepare inference summary
  summary_title <- paste0(
    "TEM-VIP (Slope) Estimate (95% CI): ",
    format(round(temvip_est, 3), digits = 3, nsmall = 3),
    " (",
    format(round(temvip_ci_lower, 3), digits = 3, nsmall = 3),
    "; ",
    format(round(temvip_ci_upper, 3), digits = 3, nsmall = 3),
    ")"
  )

  # check if the modifier is binary
  if (identical(sort(unique(modifier_vals)), c(0, 1))) {

    # plot the estimated slimple linear regression line
    ggplot2::ggplot(
      plotting_tbl,
      ggplot2::aes(x = x, y = y)
    ) +
      ggplot2::geom_point() +
      ggplot2::xlab(modifier_name) +
      ggplot2::ylab("Conditional Average Causal Effect") +
      ggplot2::ggtitle(summary_title) +
      ggplot2::theme_bw()

  } else {

    # plot the estimated slimple linear regression line
    ggplot2::ggplot(
      plotting_tbl,
      ggplot2::aes(x = x, y = y)
    ) +
      ggplot2::geom_line() +
      ggplot2::xlab(modifier_name) +
      ggplot2::ylab("Conditional Average Causal Effect") +
      ggplot2::ggtitle(summary_title) +
      ggplot2::theme_bw()
  }

}