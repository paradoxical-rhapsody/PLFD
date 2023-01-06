#' @title Detect Differential Structure
#' @description \loadmathjax
#' Identify \mjseqn{supp(M_1 - M_2)}.
#' 
#' @param x Array.
#' @param y Vector.
#' @param blockMode See [plfd()].
#' 
#' @return Logical matrix, wherein the set of `FALSE` indicates non-zero entries
#'  in \mjseqn{supp(M_1 - M_2)}.
#' 
#' @noRd
get_suppSet <- function (x, y, blockMode=NULL) {
    rDim <- NROW(x)
    cDim <- NCOL(x)
    n0 <- dim(x)[3]
    stopifnot( n0 == length(y) )
    stopifnot( y %in% 1:2 )

    if ( is.null(blockMode) )
        return( matrix(FALSE, rDim, cDim) )

    stopifnot( blockMode == "forward" )
    flag <- matrix(TRUE, rDim, cDim)
    
    logL <- cxx_prec(x[, , y==1, drop=FALSE], 
                    x[, , y==2, drop=FALSE], flag)$logLik
    # df <- (rDim+cDim)*(rDim+cDim+1) / 2.0
    ebic <- ( -2*logL +  sum(!flag) * log(n0) )

    k <- 0
    i0 <- NA
    while (any(flag)) {
        candidates <- which(flag)
        for (iC in candidates) {
            flag[iC] <- FALSE
            templogLik <- cxx_prec(x[, , y==1, drop=FALSE], 
                                    x[, , y==2, drop=FALSE], flag)$logLik
            flag[iC] <- TRUE

            if (templogLik > logL) {
                logL <- templogLik
                i0 <- iC
            }
        }
        
        flag[i0] <- FALSE
        tmpebic <- -2*logL + sum(!flag)*log(n0) # + (2*ebic.gamma*lchoose(rDim*cDim, k))
        if (tmpebic > ebic) {
            if (length(candidates) < rDim*cDim) # avoid full `TRUE`
                flag[i0] <- TRUE
            break
        }
        
        ebic <- tmpebic
    }

    return(flag)
}
