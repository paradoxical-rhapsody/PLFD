#' @title \mjeqn{T^2}{T^2} Statistic
#' @description \loadmathjax
#' 
#' @param x1 See [get_suppSet()].
#' @param x2 See [get_suppSet()].
#' 
#' @details 
#' \mjeqn{T^2 = \frac{n_1 n_2}{n} tr [(\hat{M_1} - \hat{M_2})^\top 
#'  \hat{\Sigma}^{-1} (\hat{M_1} - \hat{M_2}) \hat{Sigma}^{-1} ]}{
#'   (n1*n2/n) tr[(\hat{M1}-\hat{M2})' \hat{Psi}^{-1} (\hat{M1}-\hat{M2}) \hat{Sigma}^{-1}]
#' }
#' The total sample size should be larger than the row size and
#' column size so that the estimated row and column covariance matrices 
#' are non-singular.
#' 
#' @return \mjeqn{T^2}{T^2}-statistic.
#' @noRd
get_T2 <- function (x1, x2) {
    stopifnot(NROW(x2) == NROW(x1))
    stopifnot(NCOL(x2) == NCOL(x1))
    r0 <- NROW(x1)
    c0 <- NCOL(x1)
    n1 <- dim(x1)[3]
    n2 <- dim(x2)[3]
    n  <- n1 + n2

    M.diff <- apply(x1, 1:2, mean) - apply(x2, 1:2, mean)
    dim(M.diff) <- c(r0, c0)
    p <- cxx_prec(x1, x2, matrix(FALSE, r0, c0))

    # T2 = (n1*n2/n) tr[(M1 - M2)' Sig.inv (M1 - M2) Psi.inv]
    T2 <- (n1*n2/n) * sum(crossprod(p$PsiInv, M.diff) * tcrossprod(M.diff, p$SigInv))
    (n-r0*c0-3)/(r0*c0*(n-2)) * T2 - 1.0
}
