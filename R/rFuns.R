#' @title Scaled \eqn{T^2}-statistic
#' @param x1 Array, training data from group 1.
#' @param x2 Array, training data from group 2.
#' @return The scaled Hotelling's \eqn{T^2}-statistic.
#' @details The total samples size should be larger than the row size and the
#' column size so that the estimated row and column covariance matrices are non-singular.
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
    dim(M.diff) <- c(NROW(x1), NCOL(x1))
    p <- cxx_prec(x1, x2, matrix(FALSE, r0, c0))

    # T2 = (n1*n2/n) tr[(M1 - M2)' Sig.inv (M1 - M2) Psi.inv]
    T2 <- (n1*n2/n) * sum(crossprod(p$PsiInv, M.diff) * tcrossprod(M.diff, p$SigInv))
    (n-r0*c0-3)/(r0*c0*(n-2)) * T2 - 1.0
}

#' @title Detect the Nonzero Positions of \eqn{M_1 - M_2}
#' @description Three modes are provided for identifying \eqn{supp(M_1 - M_2)}. See details.
#' @param x1,x2 Array with common row size and column size.
#' @param blockMode \code{dense}(the default), \code{fd(forward)} or \code{bd(backward)}. See details.
#' @return Logical matrix with \code{FALSE} corresponds to the nonzero entry.
#' @details If \code{blockMode=='dense'}, the support set is assumed as all the entries. 
#' If \code{blockMode=='fd'} or \code{'bd'}, the support set is detected by forward/backward 
#' sequential likelihood method, in which BIC serves as the stoppting rule.
#' @noRd 
get_suppSet <- function (x1, x2, blockMode) {
    stopifnot(NROW(x2) == NROW(x1))
    stopifnot(NCOL(x2) == NCOL(x1))
    stopifnot(blockMode %in% c('dense', 'fd', 'bd', 'forward', 'backward'))
    if (blockMode == 'forward')  blockMode <- 'fd'
    if (blockMode == 'backward') blockMode <- 'bd'
    r0 <- NROW(x1)
    c0 <- NCOL(x1)
    n0 <- dim(x1)[3] + dim(x2)[3]

    if (blockMode == 'dense') {
        flag <- matrix(FALSE, r0, c0)
    } else {
        s0 <- switch(blockMode, 'fd'=TRUE, 'bd'=FALSE)
        flag <- matrix(s0, r0, c0)
        candidates <- 1:(r0*c0)
        
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
            df <- (r0+c0)*(r0+c0+1)/2 + sum(!flag)
            # logSize <- log_binom(r0*c0, k)
            ebicOld <- ebic
            ebic <- -2*logL + df*log(n0) #+ (2*ebic.gamma*logSize)
            
            k <- k + 1
            if (ebic > ebicOld) {
                if (k == 2) { 
                    flag <- matrix(s0, r0, c0)
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

#' @title The Parameters Involved Matrix-varite Linear Discriminant Function
#' @description Coefficient matrix and the averaged mean matrix.
#' @param x1 Array, training data from group 1.
#' @param x2 Array, training data from group 2.
#' @param blockMode The mode how the nonzero positions in \eqn{M_1 - M_2} are detected. See details.
#' @return list(M, B), in which \code{M} is the averged mean matrix and \code{B} is the coefficient 
#' matrix. See details.
#' @details Three modes are provided to \code{blockMode}, \code{dense, fd(forward), bd(backward)}. 
#' If \code{blockMode=='dense'}, the nonzero positions are assumed as all positions. If 
#' \code{blockMode=='fd'} or \code{'bd'}, the nonzero positions are detected by forward/backward 
#' sequential likelihood method, in which BIC serves as the stoppting rule.
#' The local discriminant function is \code{W(X) = tr((x[rIdx, cIdx] - M)\' B)}, where \code{M} and
#' \code{B} are to be estimated.
#' @noRd 
get_discriminantPara <- function (x1, x2, blockMode) {
    stopifnot(NROW(x1) == NROW(x2))
    stopifnot(NCOL(x1) == NCOL(x2))
    flag <- get_suppSet(x1, x2, blockMode)       
    m <- cxx_mean(x1, x2, flag)
    p <- cxx_prec(x1, x2, flag)
    B <- p[['PsiInv']] %*% (m[, , 1] - m[, , 2]) %*% p[['SigInv']]
    M <- (m[, , 1] + m[, , 2]) / 2
    return(list(M=M, B=B))
}

#' @title Split Matrix-variate into Non-overlapped Submatrices with Uniform Size
#' @param rDim,cDim Dimension of original matrix.
#' @param rSize,cSize Dimension of submatrix.
#' @return List in which each component includes the row index and the column 
#' index of one submatrix.
#' @examples 
#' size2index(30, 25, 5, 5)
#' size2index(30, 25, 4, 4)
#' @noRd 
size2blocks <- function(rDim, cDim, rSize, cSize) {
    stopifnot(rDim > rSize)
    stopifnot(cDim > cSize)
    if (rDim %% rSize) warning('rDim cannot be divided by rSize completely.')
    if (cDim %% cSize) warning('cDim cannot be divided by cSize completely.')

    k <- 0
    egg <- list()
    rBlkNum <- ceiling(rDim/rSize)
    cBlkNum <- ceiling(cDim/cSize)
    for (iR in seq(rBlkNum)) { for (iC in seq(cBlkNum)) {
        k <- k + 1
        egg[[k]] <- list(
            rIdx = (rSize*(iR-1)+1) : min(rDim, rSize*iR),
            cIdx = (cSize*(iC-1)+1) : min(cDim, cSize*iC)
        )
    }}
    return(egg)
}

#' @title Identify Sigificant Blocks
#' @description Identify the blocks whose (scaled) \eqn{T^2}-statistics are larger than a
#'  permutated threshold.
#' @param x1 Array, training data from group 1.
#' @param x2 Array, training data from group 2.
#' @param rSize Integer, the uniform row size of the submatrix-variate under consideration.
#' @param cSize Integer, the uniform column size of the submatrix-variate under consideration.
#' @param blockList List including the index of blocks.
#' @param permNum The round number of permutation that the threshold is averaged from.
#' @param alpha The upper quantile at level-\code{alpha} of the permutated statistics of 
#' blocks under consideration is used.
#' @return List with each component including the index of rows and columns of significant blocks.
#' @noRd 
get_blocks <- function (x1, x2, rSize, cSize, blockList, permNum, alpha) {
    stopifnot(NROW(x2) == NROW(x1))
    stopifnot(NCOL(x2) == NCOL(x1))
	r0 <- NROW(x1)
	c0 <- NCOL(x1)
	n1 <- dim(x1)[3]
	n2 <- dim(x2)[3]
	n  <- n1 + n2
    if (missing(blockList)) blockList <- size2blocks(r0, c0, rSize, cSize)
    
    T2 <- rep(NA_real_, length(blockList))
    for (i in seq(blockList)) {
        rIdx <- blockList[[i]][['rIdx']]
        cIdx <- blockList[[i]][['cIdx']]
        T2[i] <- get_T2(x1[rIdx, cIdx, , drop=FALSE], x2[rIdx, cIdx, , drop=FALSE])
    }

    x <- c(x1, x2)
    dim(x) <- c(r0, c0, n)
    T2Perm <- rep(NA_real_, permNum)
    for (iP in seq(permNum)) {
        temp <- rep(NA_real_, length(blockList))
        for (i in seq(blockList)) {
            rIdx <- blockList[[i]][['rIdx']]
            cIdx <- blockList[[i]][['cIdx']]
            nIdx <- sample(1:n, n1) 
            temp[i] <- get_T2(x[rIdx, cIdx, nIdx, drop=F], x[rIdx, cIdx, -nIdx, drop=F])
        }
        T2Perm[iP] <- quantile(temp, 1-alpha)
    }

    i0 <- which(T2 > mean(T2Perm))
    if (length(i0)) {
        return(blockList[i0])
    } else {
        return(blockList[which.max(T2)])
    }
}

#' @title Estimate discriminant parameters of given blocks.
#' @description Estimate all the discriminant parameters of significant blocks.
#' @param x1 Numeric Matrix of \eqn{r \times c \times n_1}.
#' @param x2 Numeric Matrix of \eqn{r \times c \times n_2}.
#' @param blockList List including the index of rows and columns of significant blocks.
#' @param blockMode The mode how the nonzero positions in \eqn{M_1 - M_2} are detected. See details.
#' @return List Each element corresponds to one block given.
#' @details Three modes are provided to \code{blockMode}, \code{dense, fd(forward), bd(backward)}. 
#' If \code{blockMode=='dense'}, the nonzero positions are assumed as all positions. If 
#' \code{blockMode=='fd'} or \code{'bd'}, the nonzero positions are detected by forward/backward 
#' sequential likelihood method, in which BIC serves as the stoppting rule.
#' @noRd 
get_paras <- function(x1, x2, blockList, blockMode) {
    for (i in seq(blockList)) {
        rIdx <- blockList[[i]][['rIdx']]
        cIdx <- blockList[[i]][['cIdx']]
        temp <- get_discriminantPara(x1[rIdx, cIdx, , drop=F], x2[rIdx, cIdx, , drop=F], blockMode)
        blockList[[i]][['B']] <- temp[['B']]
        blockList[[i]][['M']] <- temp[['M']]
    }
    return(blockList)
}

#' @title Matrix-variate Linear Discriminant Analysis
#' @description Matrix linear discriminant analysis.
#' @param x1 Array of \eqn{r \times c \times n_1}, samples from group 1.
#' @param x2 Array of \eqn{r \times c \times n_2}, samples from group 2.
#' @param rSize,cSize The common row size and column size of blocks. See details.
#' @param blockList List including the index of considered blocks. Optional if \code{rSize} and
#' \code{cSize} are provided. See details.
#' @param blockMode The mode how the nonzero positions in \eqn{M_1 - M_2} are detected. See details.
#' @param xtest (Optional)New samples to be predicted.
#' @param ytest (Optional)Vector with \eqn{1,2} entries corresponds to the labels of \code{xtest}, .
#' @param permNum The number of permutatation(default=50).
#' @param alpha The \code{alpha}-upper quantile of the permutated statistics of blocks under 
#' consideration is used(default=0).
#' @return List, \itemize{
#'  \item \code{paras} List including the parameters of significant blocks.
#'  \item \code{y} Self-predicted results for training data. It is a matrix of \eqn{(n_1+n_2)\times 2}, 
#'          the first column is the scores and the second column is the predicted labels.
#'  \item \code{mcr} The self-predicted misclassification rate for training samples.
#'  \item \code{ytest.hat} The predicted result for \code{xtest} if it is provided. It is a
#'        matrix where the first column is scores and the second column is predicted group.
#'  \item \code{mcr.test} The misclassification rate for \code{xtest} if \code{ytest} is provided.
#' }
#' @details 
#' There are two manners to specify the submatrix-variates under consideration. In the case that 
#' the matrix-variate is splitted into non-overlapped blocks that share the common row size and
#' column size, these sizes can be specified by \code{rSize} and \code{cSize}. Otherwise, the 
#' submatrix-variates can be flexibly specified by \code{blockList}, which should be a list that each 
#' component includes \code{rIdx} and \code{cIdx} corresponding to the rows index and columns 
#' index of a sumatrix-variate. See examples.
#' 
#' Three modes are provided to \code{blockMode}, \code{dense, fd(forward), bd(backward)}. 
#' If \code{blockMode=='dense'}, the nonzero positions are assumed as all positions. If 
#' \code{blockMode=='fd'} or \code{'bd'}, the nonzero positions are detected by forward/backward 
#' sequential likelihood method, in which BIC serves as the stoppting rule.
#' @examples
#' r0 <- 20
#' c0 <- 25
#' n1 <- n2 <- 100
#' ntest <- 50
#' x1 <- array(rnorm(r0*c0*n1, mean=1), dim=c(r0, c0, n1))
#' x2 <- array(rnorm(r0*c0*n2, mean=2), dim=c(r0, c0, n2))
#' xtest <- array(rnorm(r0*c0*ntest, mean=2), dim=c(r0, c0, ntest))
#' ytest <- rep(2, ntest)
#' 
#' ## Uniformly splitting
#' (result <- plfd(x1, x2, rSize=5, cSize=5, blockMode='dense', xtest=xtest, ytest=ytest))
#' 
#' ## Specify the blocks
#' blockList <- list(list(rIdx=1:5, cIdx=1:5), 
#'                   list(rIdx=6:10, cIdx=1:5), 
#'                   list(rIdx=3:9, cIdx=2:8))
#' (result <- plfd(x1, x2, blockList=blockList, blockMode='dense', xtest=xtest, ytest=ytest))
#' @importFrom stats quantile predict
#' @export 
plfd <- function(x1, x2, rSize, cSize, blockList, blockMode, 
                    xtest, ytest, permNum=50, alpha=0.0) {
    stopifnot(NROW(x1) == NROW(x2))
    stopifnot(NCOL(x1) == NCOL(x2))
    r0 <- NROW(x1)
    c0 <- NCOL(x1)
    n1 <- dim(x1)[3]
    n2 <- dim(x2)[3]
    n  <- n1 + n2
    plfd.model <- list(n1=n1, n2=n2, r0=r0, c0=c0, blockMode=blockMode, permNum=permNum, alpha=alpha)
    class(plfd.model) <- 'plfd'
    
    if (missing(blockList)) blockList <- size2blocks(r0, c0, rSize, cSize)
    plfd.model[['totalBlockNum']] <- length(blockList)

    sigBlockList <- get_blocks(x1, x2, rSize, cSize, blockList, permNum, alpha)
    paras <- get_paras(x1, x2, sigBlockList, blockMode)
    plfd.model[['paras']] <- paras
    
    y <- rbind(predict(plfd.model, x1), predict(plfd.model, x2))
    rownames(y) <- c(paste0('y1-', 1:n1), paste0('y2-', 1:n2))
    plfd.model[['y']] <- y
    plfd.model[['mcr']] <- sum(y[, 'group'] != c(rep(1, n1), rep(2, n2))) / n
    
    if (!missing(xtest)) {
        if (is.matrix(xtest)) dim(xtest) <- c(NROW(xtest), NCOL(xtest), 1)
        stopifnot(NROW(xtest) == NROW(x1))
        stopifnot(NCOL(xtest) == NCOL(x1))
        ytest.hat <- predict(plfd.model, xtest)
        plfd.model[['ytest.hat']] <- ytest.hat
    }
    if (!missing(ytest) && !missing(xtest)) {
        stopifnot(length(ytest) == dim(xtest)[3])
        plfd.model[['mcr.test']] <- sum(ytest.hat[, 'group'] != ytest) / length(ytest)
    }
    
    return(plfd.model)
}

#' @title Predict Method for \code{plfd}
#' @param object \code{plfd} object.
#' @param x The samples to be predicted.
#' @param ... Ignored currently.
#' @return Matrix with two columns, including the scores and predicted groups.
#' @examples
#' r0 <- 20
#' c0 <- 25
#' n1 <- n2 <- 100
#' x1 <- array(rnorm(r0*c0*n1, mean=1), dim=c(r0, c0, n1))
#' x2 <- array(rnorm(r0*c0*n2, mean=2), dim=c(r0, c0, n2))
#' (plfd.model <- plfd(x1, x2, rSize=5, cSize=5, blockMode='dense'))
#' 
#' ntest <- 50
#' xtest <- array(rnorm(r0*c0*ntest, mean=2), dim=c(r0, c0, ntest))
#' predict(plfd.model, xtest)
#' @export 
predict.plfd <- function(object, x, ...) {
    stopifnot(object$r0 == NROW(x))
    stopifnot(object$c0 == NCOL(x))
    if (is.matrix(x)) dim(x) <- c(NROW(x), NCOL(x), 1)
    n <- dim(x)[3]
  
    W <- rep(0.0, n)
    for (iB in 1:length(object$paras)) {
      rIdx <- object$paras[[iB]]$rIdx
      cIdx <- object$paras[[iB]]$cIdx
      M <- object$paras[[iB]]$M
      B <- object$paras[[iB]]$B
      W <- W + apply(x[rIdx, cIdx, , drop=FALSE], 3, function(.) sum((.-M) * B))
    }
  
    y <- ifelse(W>0, 1, 2)
    result <- matrix(c(W, y), n, 2)
    colnames(result) <- c('score', 'group')
    return(result)
}

#' @title Print Method for \code{plfd}
#' @param x \code{plfd} object.
#' @param ... Ignored currently.
#' @noRd 
#' @export
print.plfd <- function (x, ...) {
    cat(sprintf('Dimension of Matrix-vairate: %d x %d.\n', x$r0, x$c0))
    cat(sprintf('Consider %d block(s) in total.\n', x$totalBlockNum))
    cat(sprintf('Select %d significant block(s).\n\n', length(x$paras)))

    cat(sprintf('Training data: n1=%d, n2=%d.\n', x$n1, x$n2))
    cat(sprintf('Misclassification rate of **traning** data is %.3f.\n\n', x$mcr))

    if (!is.null(x$ytest.hat)) {
        cat('Predict testing data:\n' )
        cat(x$ytest.hat[, 'group'])
        if(!is.null(x$mcr.test)) {
            cat(sprintf('\nMisclassification rate of **testing** data is %.3f.\n', x$mcr.test))
        }
    }
}
