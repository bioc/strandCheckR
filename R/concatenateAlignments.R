#' @title Concatenate a list of Alignments into One
#' 
#' @description This method take a list of alignments and concatenate them into one.
#' 
#' @param bam a list returned by \code{scanBam} function, each element correspond to a chromosome, containing the information of strand, starting position, cigar string, and eventually flag, qname
#' @param statinfo a data frame that contains some information of the alignments
#' @param flag either the alignments contain the \code{flag} field
#' @param qname either the alignmetns contain the \code{qname} field
#' 
#' @return the concatenated alignments of the input list
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

