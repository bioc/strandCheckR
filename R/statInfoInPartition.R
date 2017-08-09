#' @title calculate the information within parition
#'
#' @param statInfo a data frame that contains some information of the alignments
#' @param winStep the winStep of sliding windows
#' 
#' @export
#'
statInfoInPartition <- function(statInfo, winStep){
  present <- which(statInfo$NbOriginalReads!=0)
  lengthInPartition <- winStep*ceiling(c(0,cumsum(as.numeric(statInfo$Length[present])))/winStep)
  nbReadsInPartition <- c(0,cumsum(statInfo$NbOriginalReads[present]))
  #update statInfo
  statInfo$FirstBaseInPartition[present] <- lengthInPartition[-length(lengthInPartition)]+1
  statInfo$LastBaseInPartition[present] <- lengthInPartition[-1]
  statInfo$FirstReadInPartition[present] <- nbReadsInPartition[-length(nbReadsInPartition)]+1
  statInfo$LastReadInPartition[present] <- nbReadsInPartition[-1]
  return(statInfo)
}
