#' @title The Parameters Involved Matrix-varite Linear Discriminant Function
#' @description Coefficient matrix and the averaged mean matrix.
#' 
#' @param x1 See [get_suppSet()].
#' @param x2 See [get_suppSet()].
#' @param blockMode See [plfd()].
#' 
#' @return `list(M, B)`, wherein `M=(M1+M2)/2` and `B = inv(Psi) * (M1 - M2) * inv(Sigma)` are 
#' the estimated parameter involved in `W(X) = tr((x - M)\' B)`.
#' @noRd
get_discriminantPara <- function (x1, x2, blockMode) {
    stopifnot(NROW(x1) == NROW(x2))
    stopifnot(NCOL(x1) == NCOL(x2))
    flag <- get_suppSet(x1, x2, blockMode)       
    m <- cxx_mean(x1, x2, flag)
    p <- cxx_prec(x1, x2, flag)
    B <- p[['PsiInv']] %*% (m[, , 1] - m[, , 2]) %*% p[['SigInv']]
    M <- (m[, , 1] + m[, , 2]) / 2
    
    list(M=M, B=B)
}
