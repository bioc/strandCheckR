#' @title Concatenate a list of Alignments into One
#' 
#' @description Concatenate a list of Alignments from multiple chromosomes into a single object
#' 
#' @details This method take a list of alignments across one or more chromosomes as output by \code{scanBam} 
#' and concatenates them into a single set of alignments which may include multiple chromosomes
#' 
#' @param bam a list returned by \code{scanBam} function, each element correspond to a chromosome, containing the information of strand, starting position, cigar string, and eventually flag, qname
#' @param statInfo a data frame that contains some key information of the alignments
#' 
#' @return the concatenated alignments of the input list
#' @export
#'
concatenateAlignments <- function(bam, statInfo){
  
  chr <- gsub("(^[^:]*):.+", "\\1", names(bam))
  nm <- names(bam[[1]])
  chk <- vapply(bam, function(x){all(names(x) %in% nm)}, logical(1))
  stopifnot(all(chk))
  stopifnot(nrow(statInfo) == length(bam))
  
  flag <- "flag" %in% nm
  qname <- "qname" %in% nm
  
  #initialize the concatenating alignments
  nbTotalReads <- sum(statInfo$NbOriginalReads)
  concatAlignments <- list("strand"=rep("*",nbTotalReads),
                           "pos"=rep(0,nbTotalReads),
                           "cigar"=rep("",nbTotalReads))
  if (flag) concatAlignments$flag <- rep(0,nbTotalReads)
  if (qname) concatAlignments$qname <- rep("",nbTotalReads)
  
  
  #concatenate the alignments
  for (i in seq_along(bam)){
    if (statInfo$NbOriginalRead[i] > 0){
      pos <- statInfo$FirstReadInPartition[i]:statInfo$LastReadInPartition[i]
      concatAlignments$strand[pos] <- as.character(bam[[i]]$strand)
      concatAlignments$pos[pos] <- bam[[i]]$pos + statInfo$FirstBaseInPartition[i]
      concatAlignments$cigar[pos] <- bam[[i]]$cigar
      if (flag) concatAlignments$flag[pos] <- bam[[i]]$flag
      if (qname) concatAlignments$qname[pos] <- bam[[i]]$qname      
    }
  }
  
  return(concatAlignments)
}

