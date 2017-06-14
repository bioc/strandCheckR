#' @title get fragments of reads from a bam file
#'
#' @export
#' @importFrom GenomicAlignments extractAlignmentRangesOnReference
#' @import S4Vectors
#'
getWinOfAlignments <- function(bam,str,win,step,limit,subset,coverage=FALSE){
  if (missing(subset)){
    index <- which(bam$strand==str)
    position <- extractAlignmentRangesOnReference(bam$cigar[index],pos=bam$pos[index]) %>% data.frame() %>% dplyr::select(-c(group_name))
    range <- IRanges(position$start,position$end,position$width)
    winrange <- getWinFromIRanges(range,win,step,limit)
    mcols(winrange) <- data.frame("alignment"=index[position$group])
  }
  else{
    index <- which(bam$strand[subset]==str)
    position <- extractAlignmentRangesOnReference(bam$cigar[subset][index],pos=bam$pos[subset][index]) %>% data.frame() %>% dplyr::select(-c(group_name))
    range <- IRanges(position$start,position$end,position$width)
    winrange <- getWinFromIRanges(range,win,step,limit)
    mcols(winrange) <- data.frame("alignment"=which(subset==TRUE)[index[position$group]])
  }
  if (coverage==TRUE){
    cov <- coverage(range)
    return(list("Win"=winrange,"Coverage"=cov))
  }
  else{
    return(list("Win"=winrange))
  }
}
