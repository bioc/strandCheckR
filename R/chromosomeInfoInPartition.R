#' @title calculate the information within parition
#' @description update the information of first/base base/read of each chromosome within the partition
#' @param chromosomeInfo a data frame that contains some information of the alignments
#' @param winStep the winStep of sliding windows
#' 
#' @export
#'
chromosomeInfoInPartition <- function(chromosomeInfo, winStep){
  
  stopifnot(isValidStatInfo(chromosomeInfo))
  
  present <- which(chromosomeInfo$NbOriginalReads!=0)
  lengthInPartition <- winStep*ceiling(c(0,cumsum(as.numeric(chromosomeInfo$Length[present])))/winStep)
  nbReadsInPartition <- c(0,cumsum(chromosomeInfo$NbOriginalReads[present]))
  #update chromosomeInfo
  chromosomeInfo$FirstBaseInPartition[present] <- lengthInPartition[-length(lengthInPartition)]+1
  chromosomeInfo$LastBaseInPartition[present] <- lengthInPartition[-1]
  chromosomeInfo$FirstReadInPartition[present] <- nbReadsInPartition[-length(nbReadsInPartition)]+1
  chromosomeInfo$LastReadInPartition[present] <- nbReadsInPartition[-1]
  return(chromosomeInfo)
}


isValidStatInfo <- function(df){
  
  reqNames <- c("Sequence", "Length", "NbOriginalReads", "FirstBaseInPartition", 
                "LastBaseInPartition","FirstReadInPartition","LastReadInPartition")
  if (all(reqNames %in% names(df))) return(TRUE)
  FALSE
  
}