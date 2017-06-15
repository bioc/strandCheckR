#' @title concatenate a list of alignments into one
#' @description concatenate a list of alignments into one
#' @param bam a list returned by \code{scanBam} function, each element correspond to a chromosome, containing the information of strand, starting position, cigar string, and eventually flag, qname
#' @param statinfo
#' @param flag
#' @param qname
#' 
#' @export
#'
concatenateAlignments <- function(bam,statinfo,flag=FALSE,qname=FALSE){
  #initialize the concatenating alignments
  nbTotalReads <- sum(statinfo$NbOriginalReads)
  concatAlignments <- list("strand"=rep("*",nbTotalReads),"pos"=rep(0,nbTotalReads),"cigar"=rep("",nbTotalReads))
  
  #concatenate the alignments
  for (i in seq_along(bam)){
    if (statinfo$NbOriginalRead[i]>0){
      concatAlignments$strand[statinfo$FirstReadInPartition[i]:statinfo$LastReadInPartition[i]] <- as.character(bam[[i]]$strand)
      concatAlignments$pos[statinfo$FirstReadInPartition[i]:statinfo$LastReadInPartition[i]] <- bam[[i]]$pos + statinfo$FirstBaseInPartition[i]
      concatAlignments$cigar[statinfo$FirstReadInPartition[i]:statinfo$LastReadInPartition[i]] <- bam[[i]]$cigar
    }
  }
  if (flag){
    concatAlignments[["flag"]] <- rep(0,nbTotalReads)
    for (i in seq_along(bam)){
      if (statinfo$NbOriginalReads[i]>0){
        concatAlignments$flag[statinfo$FirstReadInPartition[i]:statinfo$LastReadInPartition[i]] <- bam[[i]]$flag
      }
    }
  }
  if (qname){
    concatAlignments[["qname"]] <- rep("",nbTotalReads)
    for (i in seq_along(bam)){
      if (statinfo$NbOriginalReads[i]>0){
        concatAlignments$qname[statinfo$FirstReadInPartition[i]:statinfo$LastReadInPartition[i]] <- bam[[i]]$qname
      }
    }
  }
  return(concatAlignments)
}

