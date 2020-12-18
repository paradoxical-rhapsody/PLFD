#' @title Predict Method for `plfd`
#' 
#' @param object `plfd` object.
#' @param x The samples to be predicted.
#' @param ... Ignored currently.
#' 
#' @return Matrix with two columns, including the scores and predicted groups.
#' @export 
predict.plfd <- function(object, x, ...) {
    stopifnot(object$rDim == NROW(x))
    stopifnot(object$cDim == NCOL(x))
    if (is.matrix(x)) dim(x) <- c(NROW(x), NCOL(x), 1)
    n <- dim(x)[3]
  
    W <- rep(0.0, n)
    for (iB in 1:length(object$paras)) {
      rIdx <- object$paras[[iB]]$rIdx
      cIdx <- object$paras[[iB]]$cIdx
      M <- object$paras[[iB]]$M
      B <- object$paras[[iB]]$B
      W <- W + apply(x[rIdx, cIdx, , drop=FALSE], 3, function(.) sum((.-M) * B))
    }
  
    y <- ifelse(W>0, 1, 2)
    result <- matrix(c(W, y), n, 2)
    colnames(result) <- c('score', 'group')
    return(result)
}
