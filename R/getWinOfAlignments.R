#' @title get fragments of reads from a bam file
#'
#' @export
#' @importFrom GenomicAlignments extractAlignmentRangesOnReference
#' @import S4Vectors
#' @importFrom dplyr select
#'
getWinOfAlignments <- function(bam,str,win,step,limit,subset,coverage=FALSE){
  if (missing(subset)){
    index <- which(bam$strand==str)
    position <- extractAlignmentRangesOnReference(bam$cigar[index],pos=bam$pos[index]) %>% data.frame() %>% select(-c(group_name))
    maxWin <- ceiling((max(position$end)-win)/step)+1
    range <- IRanges(position$start,position$end,position$width)
    winrange <- getWinFromIRanges(range,win,step,limit,maxWin)
    mcols(winrange) <- data.frame("alignment"=index[position$group])
  }
  else{
    index <- which(bam$strand[subset]==str)
    position <- extractAlignmentRangesOnReference(bam$cigar[subset][index],pos=bam$pos[subset][index]) %>% data.frame() %>% select(-c(group_name))
    maxWin <- ceiling((max(position$end)-win)/step)+1
    range <- IRanges(position$start,position$end,position$width)
    winrange <- getWinFromIRanges(range,win,step,limit,maxWin)
    mcols(winrange) <- data.frame("alignment"=which(subset==TRUE)[index[position$group]])
  }
  if (coverage==TRUE){
    return(list("Win"=winrange,"Coverage"=coverage(range)))
  }
  else{
    return(list("Win"=winrange))
  }
}
