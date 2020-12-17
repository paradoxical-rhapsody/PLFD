#' @title Matrix-variate Linear Discriminant Analysis
#' @description Matrix linear discriminant analysis.
#' 
#' @param x1 Array of \eqn{r \times c \times n_1}, samples from group 1.
#' @param x2 Array of \eqn{r \times c \times n_2}, samples from group 2.
#' @param r0,c0 The common row size and column size of blocks. See details.
#' @param blockList List including the index of considered blocks. Optional if `r0` and
#' `c0` are provided. See details.
#' @param blockMode How the nonzero positions of `M1-M2` are detected.
#' Three modes are provided: `"dense"`(default), `"fd"/"forward"` and `"bd"/"backward"`. 
#' The `blockMode=="dense"` supposes that `M1 - M2` has no zero entry. 
#' If `blockMode=="fd"` or `"bd"`, the nonzero positions are detected by 
#' forward/backward step, where *BIC* serves as the stopping rule.
#' @param xtest (Optional) New samples to be predicted.
#' @param ytest (Optional) Vector with \eqn{1,2} entries corresponds to the labels of \code{xtest}.
#' @param permNum The number of permutation (default=50).
#' @param alpha The upper alpha` quantile of the permutation statistics (default=0). 
#' 
#' @return List, \itemize{
#'  \item \code{paras} List including the parameters of significant blocks.
#'  \item \code{y} Self-predicted results for training data. It is a matrix of \eqn{(n_1+n_2)\times 2}, 
#'          the first column is the scores and the second column is the predicted labels.
#'  \item \code{mcr} The self-predicted misclassification rate for training samples.
#'  \item \code{ytest.hat} The predicted result for \code{xtest} if it is provided. It is a
#'        matrix where the first column is scores and the second column is predicted group.
#'  \item \code{mcr.test} The misclassification rate for \code{xtest} if \code{ytest} is provided.
#' }
#' 
#' @details 
#' There are two manners to specify the blocks under consideration. In the case that 
#' the matrix-variate is partition into non-overlapping blocks that share the common row size and
#' column size, these sizes can be specified by \code{r0} and \code{c0}. Otherwise, the 
#' blocks can be flexibly specified by \code{blockList}, which should be a list that each 
#' component includes \code{rIdx} and \code{cIdx} corresponding to the rows index and columns 
#' index of a submatrix-variate. See examples.
#' 
#' @examples
#' ## Simulate the data of small dimension and sample size for saving time
#' set.seed(2020)
#' rDim <- 20
#' cDim <- 20
#' n1 <- n2 <- 50
#' ntest <- 50
#' x1 <- array(rnorm(rDim*cDim*n1, mean=0.0), dim=c(rDim, cDim, n1))
#' x2 <- array(rnorm(rDim*cDim*n2, mean=1.0), dim=c(rDim, cDim, n2))
#' xtest <- array(rnorm(rDim*cDim*ntest, mean=1.0), dim=c(rDim, cDim, ntest))
#' ytest <- rep(2, ntest)
#' 
#' ## Uniform partition
#' (result <- plfd(x1, x2, r0=5, c0=5, blockMode='dense', xtest=xtest, ytest=ytest))
#' 
#' ## Pre-specify feature blocks
#' blockList <- list(list(rIdx=1:5, cIdx=1:5), 
#'                   list(rIdx=6:10, cIdx=1:5), 
#'                   list(rIdx=3:9, cIdx=2:8))
#' (plfd.model <- plfd(x1, x2, blockList=blockList, blockMode='dense', xtest=xtest, ytest=ytest))
#' 
#' ## print
#' print(plfd.model)
#' 
#' ## Predict
#' predict(plfd.model, xtest)
#' 
#' @importFrom stats quantile predict
#' @export 
plfd <- function(x1, x2, r0, c0, blockList, blockMode='dense', xtest, ytest, permNum=50, alpha=0.0) {
    stopifnot(NROW(x1) == NROW(x2))
    stopifnot(NCOL(x1) == NCOL(x2))
    rDim <- NROW(x1)
    cDim <- NCOL(x1)
    n1 <- dim(x1)[3]
    n2 <- dim(x2)[3]
    n  <- n1 + n2
    plfd.model <- list(n1=n1, n2=n2, rDim=rDim, cDim=cDim, blockMode=blockMode, permNum=permNum, alpha=alpha)
    class(plfd.model) <- 'plfd'
    
    if (missing(blockList)) blockList <- size2blocks(rDim, cDim, r0, c0)
    plfd.model[['totalBlockNum']] <- length(blockList)

    sigBlockList <- get_feature_blocks(x1, x2, r0, c0, blockList, permNum, alpha)
    paras <- get_paras(x1, x2, sigBlockList, blockMode)
    plfd.model[['paras']] <- paras
    
    y <- rbind(predict(plfd.model, x1), predict(plfd.model, x2))
    rownames(y) <- c(paste0('y1-', 1:n1), paste0('y2-', 1:n2))
    plfd.model[['y']] <- y
    plfd.model[['mcr']] <- sum(y[, 'group'] != c(rep(1, n1), rep(2, n2))) / n
    
    if (!missing(xtest)) {
        if (is.matrix(xtest)) dim(xtest) <- c(NROW(xtest), NCOL(xtest), 1)
        stopifnot(NROW(xtest) == NROW(x1))
        stopifnot(NCOL(xtest) == NCOL(x1))
        ytest.hat <- predict(plfd.model, xtest)
        plfd.model[['ytest.hat']] <- ytest.hat
    }
    if (!missing(ytest) && !missing(xtest)) {
        stopifnot(length(ytest) == dim(xtest)[3])
        plfd.model[['mcr.test']] <- sum(ytest.hat[, 'group'] != ytest) / length(ytest)
    }
    
    return(plfd.model)
}
