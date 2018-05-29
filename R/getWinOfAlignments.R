#' @title get the window ranges of alignments
#' 
#' @description calculate the windows that contain each read fragment
#' 
#' @param readInfo a list contains the read information of one sequence
#' @param strand the considering strand
#' @param winWidth the window size
#' @param winStep the window winStep
#' @param readProp a read is considered to be included in a window if and only if at least \code{readProp} 
#' percent of it is in the window.
#' @param subset if we consider only a subset of the input reads
#' @param useCoverage either base on coverage or number of reads
#' 
#' @importFrom GenomicAlignments extractAlignmentRangesOnReference
#' @importFrom S4Vectors mcols
#' 
#' @export
getWinOfAlignments <- function(readInfo, strand, winWidth, winStep, readProp, useCoverage=FALSE, subset=NULL){
  
  if (is.null(subset)){
    index <- which(readInfo$strand==strand)
    # Read CIGAR strings to get fragment positions relative to windows
    position <- extractAlignmentRangesOnReference(readInfo$cigar[index],pos=readInfo$pos[index]) 
    position <- as.data.frame(position)
    # Get the final window number
    maxWin <- ceiling((max(position$end) - winWidth) / winStep) + 1
    range <- IRanges::IRanges(position$start, position$end, position$width)
    winrange <- getWinFromIRanges(range, winWidth, winStep, readProp, maxWin)
    mcols(winrange) <- data.frame("alignment"=index[position$group])
  }
  else{
    if (is.logical(subset) && length(subset) != length(readInfo$strand)) stop("Invalid subset vector")
    index <- which(readInfo$strand[subset]==strand)
    # Read CIGAR strings to get fragment positions relative to windows
    position <- extractAlignmentRangesOnReference(readInfo$cigar[subset][index],pos=readInfo$pos[subset][index]) 
    position <- as.data.frame(position)
    # Get the final window number
    maxWin <- ceiling((max(position$end)-winWidth) / winStep) + 1
    range <- IRanges::IRanges(position$start, position$end, position$width)
    winrange <- getWinFromIRanges(range, winWidth, winStep, readProp, maxWin)
    mcols(winrange) <- data.frame("alignment" = which(subset)[index[position$group]])
  }
  if (useCoverage==TRUE){
    return(list("Win" = winrange, "Coverage" = IRanges::coverage(range)))
  }
  else{
    return(list("Win" = winrange))
  }
}
