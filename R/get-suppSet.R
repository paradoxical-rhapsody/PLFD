#' @title Detect the Nonzero Positions of `M1 - M2`
#' @description Three modes are provided for identifying \eqn{supp(M_1 - M_2)}. See details.
#' 
#' @param x1 Array.
#' @param x2 Array with the same row and column size as `x1`.
#' @param blockMode See [plfd()].
#' 
#' @return Logical matrix with `FALSE` corresponds to the nonzero entry.
#' @noRd
get_suppSet <- function (x1, x2, blockMode) {
    stopifnot(NROW(x2) == NROW(x1))
    stopifnot(NCOL(x2) == NCOL(x1))
    stopifnot(blockMode %in% c('dense', 'fd', 'bd', 'forward', 'backward'))
    if (blockMode == 'forward')  blockMode <- 'fd'
    if (blockMode == 'backward') blockMode <- 'bd'
    rDim <- NROW(x1)
    cDim <- NCOL(x1)
    n0 <- dim(x1)[3] + dim(x2)[3]

    if (blockMode == 'dense') {
        flag <- matrix(FALSE, rDim, cDim)
    } else {
        s0 <- switch(blockMode, 'fd'=TRUE, 'bd'=FALSE)
        flag <- matrix(s0, rDim, cDim)
        candidates <- 1:(rDim*cDim)
        
        ebic <- eibcOld <- Inf
        k <- 0
        while (length(candidates)) {
            temp <- rep(NA, length(candidates))
            for (i in 1:length(candidates)) {
                flag[candidates[i]] <- !s0
                temp[i] <- cxx_prec(x1, x2, flag)$logLik
                flag[candidates[i]] <- s0
            }
            
            i0 <- candidates[switch(blockMode, 'fd'=which.max, 'bd'=which.min)(temp)]
            flag[i0] <- !s0
            logL <- cxx_prec(x1, x2, flag)$logLik
            df <- (rDim+cDim)*(rDim+cDim+1)/2 + sum(!flag)
            ebicOld <- ebic
            ebic <- -2*logL + df*log(n0) # + (2*ebic.gamma*lchoose(rDim*cDim, k))
            
            k <- k + 1
            if (ebic > ebicOld) {
                if (k == 2) { 
                    flag <- matrix(s0, rDim, cDim)
                } else { 
                    flag[i0] <- s0
                }
                break
            }
            
            candidates <- setdiff(candidates, i0)
        }
    }
    
    flag
}
