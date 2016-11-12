#' @title Filter Bam Files Based on Read Count
#' 
#' @description Filter several bamfiles based on the strand specific of all of them.
#' Strand of sliding windows are caculated based on read counts.
#' 
#' @param bamfilein the input bam files to be filterd
#' @param bamfileout the output bam files
#' @param chromosomes the list of chromosomes to be filtered
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' @param threshold the threshold upper which we keep the reads
#' @param pvalueThreshold the threshold for the p-value
#' @param minCov if a window has the max coverage least than minCov, then it will be rejected
#' @param maxCov if a window has the max coverage greater than maxCov, then it will be kept
#' @param limit the proportion that a read must be inside a window in order to be counted
#' 
#' @importFrom GenomicAlignments readGAlignments coverage 
#' @import rbamtools 
#' @import Rcpp
#' @import dplyr
#' @import magrittr
#' 
#' @return the proportion of removed reads
#' 
#' @export
filterCount <- function(bamfilein,bamfileout,chromosomes=NULL,win=1000,step=100,pvalueThreshold=0.05,minCov=0,maxCov=0,limit=0.25,threshold=0.7){
  #open a reader of the input bamfile to get chromosome names and their lengths
  reader1 <- rbamtools::bamReader(bamfilein[1],idx=TRUE) 
  refSeqs <- getRefData(reader1)
  bamClose(reader1)
  allChromosomes <- refSeqs$SN
  if (is.null(chromosomes)) chromosomes <- allChromosomes
  lenSeq <- refSeqs$LN
  #calculate the kept reads
  keep <- keepCount(bamfilein,chromosomes,allChromosomes,lenSeq,win,step,pvalueThreshold,minCov,maxCov,limit,threshold)
  for (i in seq_along(bamfilein)){#write kept reads for each bamfile
    message("Writing sample ", i)
    rbamtools::writeBam(keep[[i]],bamfilein[i],bamfileout[i],chromosomes,allChromosomes,lenSeq)
    keep[[i]] <- c(1)
  }
}
