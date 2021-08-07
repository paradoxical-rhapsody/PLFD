#' @title Print Method for `plfd`
#' 
#' @param x `plfd` object.
#' @param ... Ignored currently.
#' 
#' @export
print.plfd <- function (x, ...) {
    cat(sprintf('Dimension of Matrix-vairate: %d x %d.\n', x$rDim, x$cDim))
    cat(sprintf('Consider %d block(s) in total.\n', x$BlockNumber))
    cat(sprintf('Select %d significant block(s).\n\n', length(x$paras)))

    cat(sprintf('Training data: n1=%d, n2=%d.\n', x$n1, x$n2))
    cat(sprintf('Misclassification rate of **traning** data is %.3f.\n\n', x$mcr))

    if (!is.null(x$ytest.hat)) {
        # cat('Predict testing data:\n' )
        # cat(x$ytest.hat[, 'group'])
        if(!is.null(x$mcr.test)) {
            cat(sprintf('\nMisclassification rate of **testing** data is %.3f.\n', x$mcr.test))
        }
    }
}
