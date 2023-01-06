#' @title Estimate Discriminant Parameters for Feature Blocks.
#' 
#' @param x See [plfd()].
#' @param y See [plfd()].
#' @param blockList See [plfd()].
#' @param blockMode See [plfd()].
#' 
#' @return List. Each element corresponds to one block.
#' 
#' @noRd 
get_paras <- function(x, y, blockList, blockMode) {
    for (i in seq(blockList)) {
        rIdx <- blockList[[i]][['rIdx']]
        cIdx <- blockList[[i]][['cIdx']]
        temp <- get_discriminantPara(
                    x[rIdx, cIdx, , drop=F], y, blockMode)
        blockList[[i]][['B']] <- temp[['B']]
        blockList[[i]][['M']] <- temp[['M']]
    }
    
    return(blockList)
}
