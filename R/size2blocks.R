#' @title Split Matrix into Non-overlapped Blocks
#' 
#' @param rDim,cDim Dimension of original matrix.
#' @param r0,c0 Dimension of blocks.
#' 
#' @return List in which each component includes the row and 
#' column index set of a block.
#' 
#' @examples 
#' print( size2blocks(30, 25, 5, 5) )
#' print( size2blocks(30, 25, 4, 4) )
#' 
#' @noRd
size2blocks <- function(rDim, cDim, r0, c0){
    stopifnot(rDim >= r0)
    stopifnot(cDim >= c0)
    if (rDim %% r0) warning('rDim cannot be divided by r0 completely.')
    if (cDim %% c0) warning('cDim cannot be divided by c0 completely.')

    k <- 0
    egg <- list()
    rNum <- ceiling( rDim / r0 )
    cNum <- ceiling( cDim / c0 )
    for (iR in seq(rNum)) { for (iC in seq(cNum)) {
        k <- k + 1
        egg[[k]] <- list(
            rIdx = (r0*(iR-1)+1) : min(rDim, r0*iR),
            cIdx = (c0*(iC-1)+1) : min(cDim, c0*iC)
        )
    }}

    egg
}
