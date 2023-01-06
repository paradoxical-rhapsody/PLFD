#' @title Parameters in Discriminant Function
#' @description \loadmathjax
#' Coefficient matrix and the averaged mean matrix.
#' 
#' @param x See [get_suppSet()].
#' @param y See [get_suppSet()].
#' @param blockMode See [plfd()].
#' 
#' @return `list(M, B)`, wherein \mjseqn{M = (M_1 + M_2)/2}
#' and \mjseqn{B = Psi^{-1} (M_1 - M_2) Sigma^{-1}} are 
#' the parameters in \mjseqn{W(X) = tr((X-M)^\top B)}.
#' 
#' @noRd
get_discriminantPara <- function (x, y, blockMode) {
    stopifnot( length(dim(x)) == 3 )
    stopifnot( length(y) == dim(x)[3] )
    stopifnot( y %in% 1:2 )
    
    flag <- get_suppSet(x, y, blockMode)       
    m <- cxx_mean(x[, , y==1, drop=F], 
                x[, , y==2, drop=F], flag)
    p <- cxx_prec(x[, , y==1, drop=F], 
                x[, , y==2, drop=F], flag)
    
    list(
        M = (m$M1 + m$M2) / 2.0, 
        B = p$PsiInv %*% (m$M1 - m$M2) %*% p$SigInv
    )
}
