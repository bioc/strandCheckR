#' @title calculate the information within parition
#' @description update the information of first/base base/read of each sequence within the partition
#' @param sequenceInfo a data frame that contains some information of the alignments
#' @param winStep the winStep of sliding windows
#' 
#' @export
#'
sequenceInfoInPartition <- function(sequenceInfo, winStep){
  
  stopifnot(isValidStatInfo(sequenceInfo))
  
  present <- which(sequenceInfo$NbOriginalReads!=0)
  lengthInPartition <- winStep*ceiling(c(0,cumsum(as.numeric(sequenceInfo$Length[present])))/winStep)
  nbReadsInPartition <- c(0,cumsum(sequenceInfo$NbOriginalReads[present]))
  #update sequenceInfo
  sequenceInfo$FirstBaseInPartition[present] <- lengthInPartition[-length(lengthInPartition)]+1
  sequenceInfo$LastBaseInPartition[present] <- lengthInPartition[-1]
  sequenceInfo$FirstReadInPartition[present] <- nbReadsInPartition[-length(nbReadsInPartition)]+1
  sequenceInfo$LastReadInPartition[present] <- nbReadsInPartition[-1]
  return(sequenceInfo)
}


isValidStatInfo <- function(df){
  
  reqNames <- c("Sequence", "Length", "NbOriginalReads", "FirstBaseInPartition", 
                "LastBaseInPartition","FirstReadInPartition","LastReadInPartition")
  if (all(reqNames %in% names(df))) return(TRUE)
  FALSE
  
}