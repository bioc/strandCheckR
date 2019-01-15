#' @title get the strand information of all windows from read information
#' @description get the number of positive/negative reads of all windows from 
#' read information obtained from \code{\link{scanBam}} function
#' @param readInfo a list contains read information returned by 
#' \code{\link{scanBam}} function when read a bam file
#' @param winWidth the width of the sliding window, 1000 by default.
#' @param winStep the step length to sliding the window, 100 by default.
#' @param readProp A read is considered to be included in a window if at least
#' \code{readProp} of it is in the window. Specified as a proportion.
#' 0.5 by default.
#' @param subset an integer vector specifying the subset of reads to consider
#' @return a DataFrame object containing the number of positive/negative reads 
#' and coverage of each window sliding 
#' 
#' @seealso \code{\link{filterDNA}}, \code{\link{getStrandFromBamFile}}

getStrandFromReadInfo <- function(
    readInfo, winWidth = 1000L, winStep = 100L, readProp = 0.5, subset = NULL
    ) 
{
    winPosAlignments <- getWinIdOverlapAlignments(
        readInfo, "+", winWidth, winStep, readProp = readProp, 
        useCoverage = TRUE, subset
        )
    winNegAlignments <- getWinIdOverlapAlignments(
        readInfo, "-", winWidth, winStep, readProp = readProp, 
        useCoverage = TRUE, subset
        )
    
    # calculate strand information based on nbr of reads 
    fromNbReads <- .calculateStrandNbReads(winPosAlignments, winNegAlignments)
    
    # calculate strand information based on coverage 
    fromCoverage <- .calculateStrandCoverage(
        winPosAlignments, winNegAlignments, winWidth, winStep
        )
    
    # fill the information of the present window into the data frame to 
    # be returned
    presentWin <- which(
        as.vector(fromCoverage$CovPos > 0 | fromCoverage$CovNeg > 0) == TRUE
        )
    return(DataFrame(
        Type = Rle("",length(presentWin)), Seq = Rle("",length(presentWin)), 
        Start = presentWin, End = Rle(0,length(presentWin)), 
        NbPos = fromNbReads$NbPos[presentWin], 
        NbNeg = fromNbReads$NbNeg[presentWin], 
        CovPos = fromCoverage$CovPos[presentWin], 
        CovNeg = fromCoverage$CovNeg[presentWin], 
        MaxCoverage = fromCoverage$MaxCoverage[presentWin], 
        File = Rle("",length(presentWin))
        ))    
}
