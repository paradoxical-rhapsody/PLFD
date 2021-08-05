#' @title Estimate Discriminant Parameters for Feature Blocks.
#' 
#' @param x1 See [plfd()].
#' @param x2 See [plfd()].
#' @param blockList See [plfd()].
#' @param blockMode See [plfd()].
#' 
#' @return List. Each element corresponds to one block.
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
