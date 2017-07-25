#' @title calculate the information within parition
#'
#' @param statinfo a data frame that contains some information of the alignments
#' @param step the step of sliding windows
#' 
#' @export
#'
statInfoInPartition <- function(statinfo,step){
  present <- which(statinfo$NbOriginalReads!=0)
  lengthInPartition <- step*ceiling(c(0,cumsum(as.numeric(statinfo$Length[present])))/step)
  nbReadsInPartition <- c(0,cumsum(statinfo$NbOriginalReads[present]))
  #update statinfo
  statinfo$FirstBaseInPartition[present] <- lengthInPartition[-length(lengthInPartition)]+1
  statinfo$LastBaseInPartition[present] <- lengthInPartition[-1]
  statinfo$FirstReadInPartition[present] <- nbReadsInPartition[-length(nbReadsInPartition)]+1
  statinfo$LastReadInPartition[present] <- nbReadsInPartition[-1]
  return(statinfo)
}
