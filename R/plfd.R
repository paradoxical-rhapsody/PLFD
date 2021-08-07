#' @title A Portmanteau Local Feature Discrimination 
#'  Approach to the Classification with High-Dimensional 
#' Matrix-Variate Data
#' @description \loadmathjax
#' 
#' @param x1 Array of \mjeqn{r \times c \times n_1}{r*c*n1}, samples from group 1.
#' @param x2 Array of \mjeqn{r \times c \times n_2}{r*c*n2}, samples from group 2.
#' @param r0,c0 Row and column size of blocks. See details.
#' @param blockList List including the index set of pre-specified blocks. See details.
#' @param blockMode How the differential structure of `M1-M2` are detected.
#' The default (`blockMode=NULL`) does NOT detect the structure of feature blocks. 
#' If `blockMode="fd"`(or `"forward"`), a forward stepwise procedure is conducted to
#' detect the nonzero positions of feature blocks, wherein BIC serves as the stopping rule.
#' @param permNum Round of permutation.
#' @param alpha The upper-\mjeqn{\alpha}{alpha} quantile of the permutation statistic. 
#' 
#' @details 
#' There are two ways to specify the blocks under consideration. In the case that 
#' the matrix-variate is partition into non-overlapping blocks that share the common 
#' row size and column size, these sizes can be specified by `r0` and `c0`. Otherwise, the 
#' blocks can be flexibly specified by parameter `blockList`, which should be a list in
#' which each element includes `rIdx` and `cIdx` corresponding to the row and column index 
#' set of a block. See examples.
#' 
#' @return List, \itemize{
#'  \item `paras` List of the parameters of feature blocks.
#'  \item `y` Self-predicted results for training data. It is a matrix of 
#'      \mjeqn{(n_1+n_2)\times 2}{(n1+n2)*2}, the first column is the 
#'      scores and the second column is the predicted labels.
#'  \item `mcr` The self-predicted misclassification rate for training samples.
#'  \item `ytest.hat` The predicted result for `xtest` if it is provided. It is a
#'        matrix where the first column is scores and the second column is predicted group.
#'  \item `mcr.test` The misclassification rate for `xtest` if `ytest` is provided.
#' }
#' 
#' @examples
#' set.seed(2020)
#' rDim <- 20
#' cDim <- 20
#' 
#' n1 <- n2 <- 50
#' x1 <- array(rnorm(rDim*cDim*n1, mean=0.0), dim=c(rDim, cDim, n1))
#' x2 <- array(rnorm(rDim*cDim*n2, mean=1.0), dim=c(rDim, cDim, n2))
#'
#' ntest <- 50
#' xtest <- array(rnorm(rDim*cDim*ntest, mean=1.0), dim=c(rDim, cDim, ntest))
#' ytest <- rep(2, ntest)
#' 
#' ## Uniform partition
#' print( plfd(x1, x2, r0=5, c0=5) )
#' 
#' ## Pre-specify feature blocks
#' blockList <- list(list(rIdx=1:5, cIdx=1:5), 
#'                   list(rIdx=6:10, cIdx=1:5), 
#'                   list(rIdx=3:9, cIdx=2:8))
#' print( plfd.model <- plfd(x1, x2, blockList=blockList) )
#' 
#' ## Predict
#' predict(plfd.model, xtest)
#' 
#' @references 
#' Z. Xu, S. Luo and Z. Chen. (2021). A Portmanteau Local Feature Discrimination Approach to the 
#'  Classification with High-dimensional Matrix-variate Data. Sankhya A.
#'  doi:[10.1007/s13171-021-00255-2](https://doi.org/10.1007/s13171-021-00255-2)
#' 
#' @importFrom stats quantile predict
#' @export 
plfd <- function(x1, x2, r0, c0, blockList, blockMode=NULL, permNum=100, alpha=0.0) {
    stopifnot(NROW(x1) == NROW(x2))
    stopifnot(NCOL(x1) == NCOL(x2))
    rDim <- NROW(x1)
    cDim <- NCOL(x1)
    n1 <- dim(x1)[3]
    n2 <- dim(x2)[3]
    n  <- n1 + n2
    plfd.model <- list(n1=n1, n2=n2, rDim=rDim, cDim=cDim, 
                    blockMode=blockMode, permNum=permNum, alpha=alpha)
    class(plfd.model) <- 'plfd'
    
    if (missing(blockList)) blockList <- size2blocks(rDim, cDim, r0, c0)
    plfd.model[['BlockNumber']] <- length(blockList)

    featureBlocks <- get_feature_blocks(x1, x2, blockList, permNum, alpha)
    paras <- get_paras(x1, x2, featureBlocks, blockMode)
    plfd.model[['paras']] <- paras
    
    # y <- rbind(predict(plfd.model, x1), predict(plfd.model, x2))
    # rownames(y) <- c(paste0('y1-', 1:n1), paste0('y2-', 1:n2))
    # plfd.model[['y']] <- y
    # plfd.model[['mcr']] <- sum(y[, 'group'] != c(rep(1, n1), rep(2, n2))) / n
    
    # if (!missing(xtest)) {
    #     if (is.matrix(xtest)) dim(xtest) <- c(NROW(xtest), NCOL(xtest), 1)
    #     stopifnot(NROW(xtest) == NROW(x1))
    #     stopifnot(NCOL(xtest) == NCOL(x1))
    #     ytest.hat <- predict(plfd.model, xtest)
    #     plfd.model[['ytest.hat']] <- ytest.hat
    # }
    # if (!missing(ytest) && !missing(xtest)) {
    #     stopifnot(length(ytest) == dim(xtest)[3])
    #     plfd.model[['mcr.test']] <- sum(ytest.hat[, 'group'] != ytest) / length(ytest)
    # }
    
    return(plfd.model)
}
