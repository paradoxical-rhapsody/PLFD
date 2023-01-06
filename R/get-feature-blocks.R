#' @title Identify Feature Blocks
#' @description \loadmathjax
#' Identify the blocks whose \mjseqn{T^2}-statistic exceeds a
#'  threshold determined by permutation.
#' 
#' @param x See [plfd()].
#' @param y See [plfd()].
#' @param blockList See [plfd()].
#' @param permNum See [plfd()].
#' @param alpha See [plfd()].
#' 
#' @return List with each component including the index of rows and 
#' columns of feature blocks.
#' 
#' @noRd
get_feature_blocks <- function (x, y, blockList, permNum, alpha) {
	rDim <- NROW(x)
	cDim <- NCOL(x)
	n  <- dim(x)[3]
	n1 <- sum(y==1)
	n2 <- sum(y==2)
    stopifnot( n == length(y) )
    stopifnot( y %in% 1:2 )
    
    T2 <- rep(NA_real_, length(blockList))
    for (i in seq(blockList)) {
        rIdx <- blockList[[i]][['rIdx']]
        cIdx <- blockList[[i]][['cIdx']]
        T2[i] <- get_T2(x[rIdx, cIdx, , drop=FALSE], y)
    }

    T2Perm <- rep(NA_real_, permNum)
    for (iP in seq(permNum)) {
        temp <- rep(NA_real_, length(blockList))
        for (i in seq(blockList)) {
            rIdx <- blockList[[i]][['rIdx']]
            cIdx <- blockList[[i]][['cIdx']]
            temp[i] <- get_T2(x[rIdx, cIdx, , drop=FALSE], sample(y))
        }
        T2Perm[iP] <- quantile(temp, 1-alpha)
    }

    i0 <- union( which(T2 > mean(T2Perm)), which.max(T2) )
    return(blockList[i0])
}
