#' @title Identify Sigificant Blocks
#' @description Identify the blocks whose \eqn{T^2}-statistics are larger than a
#'  permutated thresholding.
#' 
#' @param x1 See [plfd()].
#' @param x2 See [plfd()].
#' @param r0 See [plfd()].
#' @param c0 See [plfd()].
#' @param blockList See [plfd()].
#' @param permNum See [plfd()].
#' @param alpha See [plfd()].
#' 
#' @return List with each component including the index of rows and columns of significant blocks.
#' @noRd 
get_feature_blocks <- function (x1, x2, r0, c0, blockList, permNum, alpha) {
    stopifnot(NROW(x2) == NROW(x1))
    stopifnot(NCOL(x2) == NCOL(x1))
	rDim <- NROW(x1)
	cDim <- NCOL(x1)
	n1 <- dim(x1)[3]
	n2 <- dim(x2)[3]
	n  <- n1 + n2
    if (missing(blockList)) blockList <- size2blocks(rDim, cDim, r0, c0)
    
    T2 <- rep(NA_real_, length(blockList))
    for (i in seq(blockList)) {
        rIdx <- blockList[[i]][['rIdx']]
        cIdx <- blockList[[i]][['cIdx']]
        T2[i] <- get_T2(x1[rIdx, cIdx, , drop=FALSE], x2[rIdx, cIdx, , drop=FALSE])
    }

    x <- c(x1, x2)
    dim(x) <- c(rDim, cDim, n)
    T2Perm <- rep(NA_real_, permNum)
    for (iP in seq(permNum)) {
        temp <- rep(NA_real_, length(blockList))
        for (i in seq(blockList)) {
            rIdx <- blockList[[i]][['rIdx']]
            cIdx <- blockList[[i]][['cIdx']]
            nIdx <- sample(1:n, n1)
            temp[i] <- get_T2(x[rIdx, cIdx, nIdx, drop=F], x[rIdx, cIdx, -nIdx, drop=F])
        }
        T2Perm[iP] <- quantile(temp, 1-alpha)
    }

    i0 <- which(T2 > mean(T2Perm))
    if (length(i0)) {
        return(blockList[i0])
    } else {
        return(blockList[which.max(T2)])
    }
}