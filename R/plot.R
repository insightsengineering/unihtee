#' Plot a \code{unihtee} Object
#'
#' @description \code{plot.unihtee()} produces a plot depicting the estimated
#'   simple linear regression line for the specified \code{modifier} parameter
#'   estimated by \code{\link{unihtee}()}.
#'
#' @param x A \code{unihtee} class object output by \code{\link{unihtee}()}.
#' @param ... Ignored.
#' @param modifier A \code{character} specifying the treatment effect modifier
#'   variable importance parameter estimate to plot.
#' @param print_interpretation A \cpde{flag} indicating whether to print the
#'   interpretation of the TEM-VIP inference plot. Defaults to \code{TRUE}.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_hline
#'   scale_linetype_manual labs xlab ylab ggtitle theme_bw theme
#'
#' @returns A \code{\link[ggplot2]{ggplot2}} object.
#'
#' @exportS3Method unihtee::plot
plot.unihtee <- function(
  x,
  ...,
  modifier_name,
  print_interpretation = TRUE
) {

  # extract the ACE estimate
  ace_estimate <- x$ace_estimate

  # extract the TEM-VIP inference
  temvip_inference <- x$temvip_inference_tbl[
    modifier == modifier_name, ]
  temvip_est <- temvip_inference$estimate
  temvip_ci_lower <- temvip_inference$ci_lower
  temvip_ci_upper <- temvip_inference$ci_upper

  # extract the modifier values from the data
  modifier_vals <- x$data[[modifier_name]]

  # compute the mean modifier value
  modifier_mean <- mean(modifier_vals)

  # create the data.frame for plotting
  plotting_tbl <- data.frame(
    x = modifier_vals,
    y = ace_estimate + temvip_est * (modifier_vals - modifier_mean),
    ace_estimate = ace_estimate
  )

  # prepare TEM-VIP estimate and nominal 95% CI summary for the plot title
  n_dig <- 3
  plot_title <- paste0(
    modifier_name, " TEM-VIP Inference [Estimate (Nominal 95% CI)]: ",
    format(round(temvip_est, n_dig), digits = n_dig, nsmall = n_dig),
    " (",
    format(round(temvip_ci_lower, n_dig), digits = n_dig, nsmall = n_dig),
    "; ",
    format(round(temvip_ci_upper, n_dig), digits = n_dig, nsmall = n_dig),
    ")"
  )

  # check if the modifier is binary
  if (identical(sort(unique(modifier_vals)), c(0, 1))) {

    # print interpretation
    if (print_interpretation) {
      cat(
        paste(
          "The TEM-VIP estimate corresponds to the vertical distance between",
          "both points."
        )
      )
    }

    # plot the estimated intercepts for each subgroupt
    ggplot2::ggplot(
      plotting_tbl,
      ggplot2::aes(x = factor(x, levels = c(0, 1)), y = y)
    ) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(
        ggplot2::aes(
          yintercept = ace_estimate,
          linetype = "Average Causal Effect Estimate"
        ),
        colour = "red"
      ) +
      ggplot2::scale_linetype_manual(values = 2) +
      ggplot2::labs(linetype = NULL) +
      ggplot2::xlab(modifier_name) +
      ggplot2::ylab("Conditional Average Causal Effect") +
      ggplot2::ggtitle(plot_title) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position="bottom")

  } else {

    # print interpretation
    if (print_interpretation) {
      cat(
        paste(
          "The TEM-VIP estimate corresponds to the slope of the black line."
        )
      )
    }

    # plot the estimated slimple linear regression line
    ggplot2::ggplot(
      plotting_tbl,
      ggplot2::aes(x = x, y = y)
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_hline(
        ggplot2::aes(
          yintercept = ace_estimate,
          linetype = "Average Causal Effect Estimate"
        ),
        colour = "red"
      ) +
      ggplot2::scale_linetype_manual(values = 2) +
      ggplot2::labs(linetype = NULL) +
      ggplot2::xlab(modifier_name) +
      ggplot2::ylab("Conditional Average Causal Effect") +
      ggplot2::ggtitle(plot_title) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position="bottom")
  }

}