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
filterMultiCount <- function(bamfilein,bamfileout,chromosomes=NULL,win=1000,step=100,threshold,pvalueThreshold=0.05,minR=0,limit=0.25){
  library(GenomicAlignments)
  library(rbamtools)
  library(Rcpp)
  library(dplyr)
  library(magrittr)
  reader1 <- bamReader(bamfilein[1],idx=TRUE) #open a reader of the input bamfile to extract read afterward
  refSeqs <- getRefData(reader1)
  remove(reader1)
  allChromosomes <- refSeqs$SN
  lenSeq <- refSeqs$LN
  keep <- keepMultiCount(bamfilein,bamfileout,chromosomes,allChromosomes,lenSeq,win,step,threshold,pvalueThreshold,minR,limit)
  writeBam(keep,bamfilein,bamfileout,chromosomes,allChromosomes,lenSeq)
}
