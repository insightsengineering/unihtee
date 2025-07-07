#' Print a \code{unihtee} Object
#'
#' @param x A \code{unihtee} object output by \code{\link{unihtee}()}
#' @param ... Ignored.
#'
#' @exportS3Method unihtee::print
print.unihtee <- function(x, ...) {
  print(x$temvip_inference_tbl, ...)
}
