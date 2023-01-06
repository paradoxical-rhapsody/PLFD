#' @title PLFD
#' @description \loadmathjax
#' A portmanteau local feature discrimination (PLFD) approach to the classification with
#' high-dimensional matrix-variate data.
#' 
#' @param x Array of \mjseqn{r \times c \times n}.
#' @param y Vector of length-\mjseqn{n} with values 1 or 2.
#' @param r0,c0 Row and column size of blocks. See details.
#' @param blockList List including the index set of pre-specified blocks. See details.
#' @param blockMode How the differential structure of \mjseqn{M_1 - M_2} are 
#' detected. The default (`blockMode=NULL`) does NOT detect the structure of feature 
#' blocks. If `blockMode="fd"`(or `"forward"`), a forward stepwise procedure is 
#' conducted to detect the nonzero positions of feature blocks, wherein BIC serves 
#' as the stopping rule.
#' @param permNum Rounds of permutation.
#' @param alpha The **upper**-\mjseqn{\alpha} quantile of the permutation statistic. 
#' 
#' @details 
#' There are two ways to specify the blocks under consideration. In the case that 
#' the matrix-variate is partition into non-overlapping blocks that share the common 
#' row size and column size, these sizes can be specified by `r0` and `c0`. Otherwise, the 
#' blocks can be flexibly specified by parameter `blockList`, which should be a list in
#' which each element includes `rIdx` and `cIdx` corresponding to the row and column index 
#' set of a block. See examples.
#' 
#' @return List.
#'  * `n1`, `n2`, `rDim`, `cDim`, `blockMode`, `permNum`, `alpha`;
#'  * `blockNumber`: the number of identified feature blocks.
#'  * `paras`: `list(list(rIdx, cIdx, B, M), ...)`, list of the information of 
#'  feature blocks.
#' 
#' @examples
#' set.seed(2023)
#' rDim <- 20
#' cDim <- 20
#' 
#' n <- 100
#' y <- sample(1:2, n, TRUE, c(0.5, 0.5))
#' x <- array(rnorm(rDim*cDim*n), dim=c(rDim, cDim, n))
#' x[, , y==2] <- (x[, , y==2] + 1.0)
#'
#' ntest <- 200
#' ytest <- sample(1:2, ntest, TRUE, c(0.5, 0.5))
#' xtest <- array(rnorm(rDim*cDim*ntest), dim=c(rDim, cDim, ntest))
#' xtest[, , ytest==2] <- (xtest[, , ytest==2] + 1.0)
#' 
#' ## Uniform partition
#' print( plfd(x, y, r0=5, c0=5) )
#' 
#' ## Pre-specify feature blocks
#' blockList <- list(list(rIdx=1:5, cIdx=1:5), 
#'                   list(rIdx=6:10, cIdx=1:5), 
#'                   list(rIdx=3:9, cIdx=2:8))
#' print( plfd.model <- plfd(x, y, blockList=blockList) )
#' 
#' ## Predict
#' predict(plfd.model, xtest, ytest)
#' 
#' @references 
#' Xu Z., Luo S. and Chen Z. (2021). A Portmanteau Local Feature Discrimination
#'  Approach to the Classification with High-dimensional Matrix-variate Data. Sankhya A.
#'  \doi{10.1007/s13171-021-00255-2}
#' 
#' @export 
plfd <- function(x, y, r0, c0, blockList, blockMode=NULL, permNum=100, alpha=0.0) {
    rDim <- NROW(x)
    cDim <- NCOL(x)
    n  <- dim(x)[3]
    n1 <- sum( y == 1 )
    n2 <- sum( y == 2 )
    stopifnot( y %in% 1:2 )
    stopifnot( n == length(y) )

    plfd.model <- list(n1=n1, n2=n2, rDim=rDim, cDim=cDim, 
                    blockMode=blockMode, permNum=permNum, alpha=alpha)
    class(plfd.model) <- 'plfd'
    
    if (missing(blockList)) 
        blockList <- size2blocks(rDim, cDim, r0, c0)
    
    plfd.model[['BlockNumber']] <- length(blockList)
    featureBlocks <- get_feature_blocks(x, y, blockList, permNum, alpha)
    plfd.model[['paras']] <- get_paras(x, y, featureBlocks, blockMode)
    
    return(plfd.model)
}
