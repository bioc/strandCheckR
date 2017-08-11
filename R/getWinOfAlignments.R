#' @title get the window ranges of alignments
#' @description calculate the windows that contain each read fragment
#' @param bam a list contains the read information of one chromosome
#' @param str the considering strand
#' @param winSize the window size
#' @param winStep the window winStep
#' @param limit a read is considered to be included in a window if and only if at least \code{limit} percent of it is in the window.
#' @param subset if we consider only a subset of the input reads
#' @param useCoverage either base on coverage or number of reads
#' @export
#' @importFrom GenomicAlignments extractAlignmentRangesOnReference
#' @import S4Vectors
#'
getWinOfAlignments <- function(bam, str, winSize, winStep, limit, useCoverage=FALSE, subset){
  if (missing(subset)){
    index <- which(bam$strand==str)
    position <- extractAlignmentRangesOnReference(bam$cigar[index],pos=bam$pos[index]) %>% data.frame() %>% dplyr::select(-c(group_name))
    maxWin <- ceiling((max(position$end)-winSize)/winStep)+1
    range <- IRanges(position$start,position$end,position$width)
    winrange <- getWinFromIRanges(range,winSize,winStep,limit,maxWin)
    mcols(winrange) <- data.frame("alignment"=index[position$group])
  }
  else{
    index <- which(bam$strand[subset]==str)
    position <- extractAlignmentRangesOnReference(bam$cigar[subset][index],pos=bam$pos[subset][index]) %>% data.frame() %>% dplyr::select(-c(group_name))
    maxWin <- ceiling((max(position$end)-winSize)/winStep)+1
    range <- IRanges(position$start,position$end,position$width)
    winrange <- getWinFromIRanges(range,winSize,winStep,limit,maxWin)
    mcols(winrange) <- data.frame("alignment"=which(subset==TRUE)[index[position$group]])
  }
  if (useCoverage==TRUE){
    return(list("Win"=winrange,"Coverage"=coverage(range)))
  }
  else{
    return(list("Win"=winrange))
  }
}
