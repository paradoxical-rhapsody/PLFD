#' @title Print Method for `plfd`
#' 
#' @param x `plfd` object.
#' @param ... Ignored currently.
#' 
#' @export
print.plfd <- function (x, ...) {
    message(sprintf('Dimension of Matrix-variate: %d x %d.', x$rDim, x$cDim))
    message(sprintf('Training data: n1=%d, n2=%d.', x$n1, x$n2))
    message(sprintf('Number of feature block(s): %d.', length(x$paras)))
}
