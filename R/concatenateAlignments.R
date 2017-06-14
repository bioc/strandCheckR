#' @title concatenate a list of alignments into one
#'
#' @export
#'
concatenateAlignments <- function(bam,nbOriginalReadsInPart,lengthSeqInPart,nbAlignments,flag=FALSE,qname=FALSE){
  bamPart <- list("strand"=rep("*",nbAlignments),"pos"=rep(0,nbAlignments),"cigar"=rep("",nbAlignments))
  if (flag) {
    bamPart[["flag"]] <- rep(0,nbAlignments)
  }
  if (qname) {
    bamPart[["qname"]] <- rep("",nbAlignments)
  }
  for (i in 1:length(bam)){
    bamPart$strand[(nbOriginalReadsInPart[i]+1):nbOriginalReadsInPart[i+1]] <- as.character(bam[[i]]$strand)
    bamPart$pos[(nbOriginalReadsInPart[i]+1):nbOriginalReadsInPart[i+1]] <- bam[[i]]$pos + lengthSeqInPart[i]
    bamPart$cigar[(nbOriginalReadsInPart[i]+1):nbOriginalReadsInPart[i+1]] <- bam[[i]]$cigar
    if (flag){
      bamPart$flag[(nbOriginalReadsInPart[i]+1):nbOriginalReadsInPart[i+1]] <- bam[[i]]$flag
    }
    if (qname){
      bamPart$qname[(nbOriginalReadsInPart[i]+1):nbOriginalReadsInPart[i+1]] <- bam[[i]]$qname
    }
  }
  return(bamPart)
}
