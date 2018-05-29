#' @title Concatenate a list of Alignments into One
#' 
#' @description Concatenate a list of Alignments from multiple sequences into a single object
#' 
#' @details This method take a list of alignments across one or more sequences as output by \code{scanBam} 
#' and concatenates them into a single set of alignments which may include multiple sequences
#' 
#' @param bam a list returned by \code{scanBam} function, each element correspond to a sequence, containing the information of strand, starting position, cigar string, and eventually flag, qname
#' @param sequenceInfo a data frame that contains some key information of the alignments
#' 
#' @return the concatenated alignments of the input list
#' @export
#'
concatenateAlignments <- function(bam, sequenceInfo){
  
  nm <- names(bam[[1]])
  chk <- vapply(bam, function(x){all(names(x) %in% nm)}, logical(1))
  stopifnot(all(chk))
  stopifnot(nrow(sequenceInfo) == length(bam))
  
  #initialize the concatenating alignments
  nbTotalReads <- sum(sequenceInfo$NbOriginalReads)
  concatAlignments <- vector("list",length(nm))
  names(concatAlignments) <- nm
  for (name in nm){
    if (is.factor(bam[[1]][[name]])){
      concatAlignments[[name]] <- factor(rep("*",nbTotalReads), levels = levels(bam[[1]][[name]]))
    } else{
      if (typeof(bam[[1]][[name]]) == "integer"){
        concatAlignments[[name]] <- rep(0,nbTotalReads)
      } 
      if (typeof(bam[[1]][[name]]) == "character"){
        concatAlignments[[name]] <- rep("",nbTotalReads)
      }  
    }
  }
  #concatenate the alignments
  for (i in seq_along(bam)){
    if (sequenceInfo$NbOriginalRead[i] > 0){
      range <- sequenceInfo$FirstReadInPartition[i]:sequenceInfo$LastReadInPartition[i]
      concatAlignments$pos[range] <- bam[[i]]$pos + sequenceInfo$FirstBaseInPartition[i]
      for (name in nm[nm!="pos"]){
        concatAlignments[[name]][range] <- bam[[i]][[name]]
      }
    }
  }
  return(concatAlignments)
}

