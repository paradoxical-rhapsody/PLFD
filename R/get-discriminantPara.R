#' @title Parameters in Discriminant Function
#' @description \loadmathjax
#' Coefficient matrix and the averaged mean matrix.
#' 
#' @param x1 See [get_suppSet()].
#' @param x2 See [get_suppSet()].
#' @param blockMode See [plfd()].
#' 
#' @return `list(M, B)`, wherein \mjseqn{M = (M_1 + M_2)/2}
#' and \mjseqn{B = Psi^{-1} (M_1 - M_2) Sigma^{-1}} are 
#' the parameters in \mjseqn{W(X) = tr((X-M)^\top B)}.
#' 
#' @noRd
get_discriminantPara <- function (x1, x2, blockMode) {
    stopifnot(NROW(x1) == NROW(x2))
    stopifnot(NCOL(x1) == NCOL(x2))
    
    flag <- get_suppSet(x1, x2, blockMode)       
    m <- cxx_mean(x1, x2, flag)
    p <- cxx_prec(x1, x2, flag)
    
    list(
        M = (m$M1 + m$M2) / 2.0, 
        B = p$PsiInv %*% (m$M1 - m$M2) %*% p$SigInv
    )
}
