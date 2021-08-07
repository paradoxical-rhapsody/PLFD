#' @title Print Method for `plfd`
#' 
#' @param x `plfd` object.
#' @param ... Ignored currently.
#' 
#' @export
print.plfd <- function (x, ...) {
    cat(sprintf('Dimension of Matrix-vairate: %d x %d.\n', x$rDim, x$cDim))
    cat(sprintf('Training data: n1=%d, n2=%d.\n', x$n1, x$n2))
    cat(sprintf('Number of feature block(s): %d.\n\n', length(x$paras)))
}
