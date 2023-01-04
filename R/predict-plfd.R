#' @title Predict Method for `plfd`
#' 
#' @param object `plfd` object.
#' @param x Array, matrix-variate data to be predicted.
#' @param y Vector (optional), Labels of `x` with value `1` or `2`.
#' @param ... Ignored currently.
#' 
#' @return `list(W, y.hat, mcr)`, \itemize{
#'  \item `W`: discriminant scores;
#'  \item `y.hat`: predicted labels;
#'  \item `mcr`: misclassification rate if parameter `y` is available.
#' }
#' @export 
predict.plfd <- function(object, x, y, ...) {
    stopifnot( object$rDim == NROW(x) )
    stopifnot( object$cDim == NCOL(x) )

    if (is.matrix(x)) dim(x) <- c(NROW(x), NCOL(x), 1)
    n <- dim(x)[3]
  
    W <- rep(0.0, n)
    for (iB in seq(object$paras)) {
      rIdx <- object$paras[[iB]]$rIdx
      cIdx <- object$paras[[iB]]$cIdx
      M <- object$paras[[iB]]$M
      B <- object$paras[[iB]]$B
      W <- W + apply(x[rIdx, cIdx, , drop=FALSE], 3, function(.) sum((.-M) * B))
    }

    y.hat <- ifelse(W>0, 1, 2)
    result <- list(score=W, y.hat=y.hat)

    if (!missing(y)) {
      stopifnot( n == length(y) ) 
      stopifnot( all(unique(y) %in%  1:2) )
      result["mcr"] = sum(y.hat != y) / n
    }

    return(result)
}
