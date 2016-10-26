#' filter several bamfiles based on the strand specific of all of them. Strand of sliding windows are caculated based on read counts.
#' @param bamfilein the input bam files to be filterd
#' @param bamfileout the output bam files
#' @param chromosomes the list of chromosomes to be filtered
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' @param threshold the threshold upper which we keep the reads
#' @param pvalueThreshold the threshold for the p-value
#' @param limit the proportion of a read that it should not exceed to be considered to be in a window
#' @export
filterCount <- function(bamfilein,bamfileout,chromosomes=NULL,win=1000,step=100,pvalueThreshold=0.05,minR=0,maxR=0,limit=0.25,threshold){
  #open a reader of the input bamfile to get chromosomes list and their lengths
  library(rbamtools)
  reader1 <- bamReader(bamfilein[1],idx=TRUE) 
  refSeqs <- getRefData(reader1)
  bamClose(reader1)
  allChromosomes <- refSeqs$SN
  if (is.null(chromosomes)) chromosomes<-allChromosomes
  lenSeq <- refSeqs$LN
  #calculate the kept reads
  keep <- keepCount(bamfilein,bamfileout,chromosomes,allChromosomes,lenSeq,win,step,pvalueThreshold,minR,maxR,limit,threshold)
  for (i in c(1:length(bamfilein))){
    message("Writing sample ", i)
    writeBam(keep[[i]],bamfilein[i],bamfileout[i],chromosomes,allChromosomes,lenSeq)
    keep[[i]] <- c(1)
  }
}
