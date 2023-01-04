#' @title Detect Differential Structure
#' @description \loadmathjax
#' Identify \mjeqn{supp(M_1 - M_2)}{supp(M1 - M2)}.
#' 
#' @param x1 Array.
#' @param x2 Array.
#' @param blockMode See [plfd()].
#' 
#' @return Logical matrix, wherein the set of `FALSE` corresponds to 
#'  \mjeqn{supp(M_1 - M_2)}{supp(M1 - M2)}.
#' @noRd
get_suppSet <- function (x1, x2, blockMode=NULL) {
    rDim <- NROW(x1)
    cDim <- NCOL(x1)
    n0 <- dim(x1)[3] + dim(x2)[3]

    stopifnot( NROW(x2) == rDim )
    stopifnot( NCOL(x2) == cDim )

    if ( is.null(blockMode) ) {
        return( matrix(FALSE, rDim, cDim) )
    } else {
        stopifnot( blockMode %in% c("fd", "forward") )
    }

    flag <- matrix(TRUE, rDim, cDim)
    
    logL <- cxx_prec(x1, x2, flag)$logLik
    df <- (rDim+cDim)*(rDim+cDim+1)/2
    ebic <- ( -2*logL +  df * log(n0))

    k <- 0
    candidates <- seq(rDim*cDim)
    while (length(candidates)) {
        temp <- rep(NA, length(candidates))
        for (i in 1:length(candidates)) {
            flag[candidates[i]] <- FALSE
            temp[i] <- cxx_prec(x1, x2, flag)$logLik
            flag[candidates[i]] <- TRUE
        }
        
        i0 <- candidates[which.max(temp)]
        flag[i0] <- FALSE
        logL <- cxx_prec(x1, x2, flag)$logLik
        df <- df + 1
        tmp <- -2*logL + df*log(n0) # + (2*ebic.gamma*lchoose(rDim*cDim, k))
        
        k <- k + 1
        if (tmp > ebic) {
            # flag[i0] <- TRUE
            break
        }
        
        candidates <- setdiff(candidates, i0)
        ebic <- tmp
    }

    flag
}
