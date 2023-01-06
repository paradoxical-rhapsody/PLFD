#' @title \mjseqn{T^2}-Statistic
#' @description \loadmathjax
#' 
#' @param x See [get_suppSet()].
#' @param y See [get_suppSet()].
#' 
#' @details 
#' \mjsdeqn{T^2 = \frac{n_1 n_2}{n} tr [(\hat{M_1} - \hat{M_2})^\top 
#'  \hat{\Sigma}^{-1} (\hat{M_1} - \hat{M_2}) \hat{Sigma}^{-1} ]}
#' 
#' The total sample size should be larger than the row size and
#' column size so that the estimated row and column covariance matrices 
#' are non-singular.
#' 
#' @return \mjseqn{T^2}-statistic.
#' 
#' @noRd
get_T2 <- function (x, y) {
    r0 <- NROW(x)
    c0 <- NCOL(x)
    n  <- dim(x)[3]
    n1 <- sum(y==1)
    n2 <- sum(y==2)
    stopifnot( n == length(y) )

    M1 <- apply(x[, , y==1, drop=FALSE], 1:2, mean)
    M2 <- apply(x[, , y==2, drop=FALSE], 1:2, mean)
    M.diff <- M1 - M2
    stopifnot( dim(M.diff) == c(r0, c0) )
    p <- cxx_prec(x[, , y==1, drop=FALSE], 
                x[, , y==2, drop=FALSE], 
                matrix(FALSE, r0, c0))

    # T2 = (n1*n2/n) tr[PsiInv (M1 - M2) SigInv (M1 - M2)']
    T2 <- (n1*n2/n) * sum( (p$PsiInv %*% M.diff) * (M.diff %*% p$SigInv) )

    return( (n-r0*c0-3)/(r0*c0*(n-2)) * T2 - 1.0 )
}
