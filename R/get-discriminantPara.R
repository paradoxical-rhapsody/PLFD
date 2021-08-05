#' @title Parameters in Discriminant Function
#' @description \loadmathjax
#' Coefficient matrix and the averaged mean matrix.
#' 
#' @param x1 See [get_suppSet()].
#' @param x2 See [get_suppSet()].
#' @param blockMode See [plfd()].
#' 
#' @return `list(M, B)`, wherein \mjeqn{M = (M_1 + M_2)/2}{`M=(M1+M2)/2`}
#' and \mjeqn{B = Psi^{-1} (M_1 - M_2) Sigma^{-1}}{`B = inv(Psi) * (M1 - M2) * inv(Sigma)`} are 
#' the parameters in \mjeqn{W(X) = tr((X-M)^\top B)}{`W(X) = tr((x - M)\' B)`}.
#' @noRd
get_discriminantPara <- function (x1, x2, blockMode) {
    stopifnot(NROW(x1) == NROW(x2))
    stopifnot(NCOL(x1) == NCOL(x2))
    
    flag <- get_suppSet(x1, x2, blockMode)       
    m <- cxx_mean(x1, x2, flag)
    p <- cxx_prec(x1, x2, flag)
    B <- p[['PsiInv']] %*% (m[['M1']] - m[['M2']]) %*% p[['SigInv']]
    M <- (m[['M1']] + m[['M2']]) / 2
    
    list(M=M, B=B)
}
