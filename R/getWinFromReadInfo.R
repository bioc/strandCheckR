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
#' @seealso \code{\link{filterDNA}}, \code{\link{getWinFromBamFile}}
#' @export
#' @examples
#' library(Rsamtools)
#' file <- system.file("extdata","s1.sorted.bam",package = "strandCheckR") 
#' readInfo <- scanBam(file)[[1]]
#' win <- getWinFromReadInfo(readInfo)

getWinFromReadInfo <- function(readInfo, winWidth = 1000, winStep = 100, 
                                readProp = 0.5, subset=NULL){
    
    winPositiveAlignments <- getWinOfAlignments(readInfo, "+", winWidth, 
        winStep, readProp = readProp, useCoverage = TRUE, subset)
    winNegativeAlignments <- getWinOfAlignments(readInfo, "-", winWidth, 
        winStep, readProp = readProp, useCoverage = TRUE, subset)
    
    ######################################################
    # calculate strand information based on nbr of reads #
    ######################################################
    fromNbReads <- calculateStrandNbReads(winPositiveAlignments,
                                        winNegativeAlignments)
    
    ##################################################
    # calculate strand information based on coverage #
    ##################################################
    fromCoverage <- calculateStrandCoverage(winPositiveAlignments,
        winNegativeAlignments,winWidth,winStep)
    
    # fill the information of the present window into the data 
    # frame to be returned
    presentWin <- which(as.vector(fromCoverage$CovPositive>0 |
                                fromCoverage$CovNegative>0)==TRUE)
    return(DataFrame(Type = "", Seq = "", 
                    Start = presentWin, End = 0,
                    NbPositive = fromNbReads$NbPositive[presentWin], 
                    NbNegative = fromNbReads$NbNegative[presentWin],
                    CovPositive = fromCoverage$CovPositive[presentWin], 
                    CovNegative = fromCoverage$CovNegative[presentWin],
                    MaxCoverage = fromCoverage$MaxCoverage[presentWin],
                    File = ""))
}
