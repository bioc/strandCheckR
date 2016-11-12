#' @title Filter Using Merged Files
#' 
#' @description Filter several bamfiles based on the strand specific of the merged files. 
#' The strand of each sliding window is caculated based on coverage.
#' 
#' @details None
#' 
#' @param bamfilein the input bam files to be filterd
#' @param bamfileout the output bam files
#' @param chromosomes the list of chromosomes to be filtered
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' @param threshold the threshold upper which we keep the reads
#' @param pvalueThreshold the threshold for the p-value
#' @param minR if a window has the max coverage least than minR, then it will be rejected
#' @param maxR if a window has the max coverage greater than maxR, then it will be kept
#'
#' @importFrom GenomicAlignments readGAlignments coverage 
#' @import rbamtools
#' @import Rcpp
#' @import dplyr
#' @import magrittr 
#' 
#' @return the proportion of removed reads in every sample
#' 
filterMultiCov <- function(bamfilein,bamfileout,chromosomes=NULL,win=1000,step=100,pvalueThreshold=0.05,minR=0,maxR=0,threshold=0.7){

  logitThreshold <- binomial()$linkfun(threshold)
  alignment <- readGAlignments(bamfilein[1]) # read the first input alignment
  covPos <- alignment[strand(alignment)=="+"] %>% coverage()
  covNeg <- alignment[strand(alignment)=="-"] %>% coverage()
  idSeq <- levels(seqnames(alignment)) #get the chromosome list
  remove(alignment)
  if (is.null(chromosomes)) chromosomes <- idSeq
  lenSeq<-sapply(1:length(covPos),function(i) length(covPos[[i]])) #get the length of each chromosome as the number of bases
  if (length(bamfilein)>1){
    for (i in c(2:length(bamfilein))){
      alignment <- readGAlignments(bamfilein[i])
      covPos <-covPos +alignment[strand(alignment)=="+"] %>% coverage()
      covNeg <-covNeg+alignment[strand(alignment)=="-"] %>% coverage()
      remove(alignment)
    }
  }
  keepWinPos <- list()
  keepWinNeg <- list()
  for (chr in chromosomes){
    chromosomeIndex <- which(idSeq==chr)
    len <- lenSeq[chromosomeIndex]
    #compute the normalized value of each positive/negative window to be tested by pnorm
    windows <- computeWinCov(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),len,win,step,minR,maxR,logitThreshold)
    windows$Plus <- mutate(windows$Plus,"pvalue"=pnorm(value,lower.tail = FALSE))
    windows$Minus <- mutate(windows$Minus,"pvalue"=pnorm(value,lower.tail = FALSE))
    keepWinPos[[chromosomeIndex]] <- filter(windows$Plus, pvalue <= pvalueThreshold)$win
    keepWinNeg[[chromosomeIndex]] <- filter(windows$Minus, pvalue <= pvalueThreshold)$win
    remove(windows)
  }
  remove(covPos)
  remove(covNeg)
  levelConta <- c()
  #write the kept reads to output files
  # seq_along(bamfilein)
  for (i in seq_along(bamfilein)){
    message("File ",bamfilein[i])
    alignment <- readGAlignments(bamfilein[i])
    reader <- bamReader(bamfilein[i],idx=TRUE)
    header <- getHeader(reader)
    writer <- bamWriter(header,bamfileout[i])
    nbOReads <- 0
    nbKReads <- 0
    for (chr in chromosomes){
      message("Chromosome ",chr)
      chromosomeIndex <- which(idSeq==chr)
      end <- lenSeq[chromosomeIndex]
      alignmentInChr <- alignment[seqnames(alignment)==chr]#get the reads in the considering chromosome
      alignment <- alignment[seqnames(alignment)!=chr]#reduce the size of alignment (for memory purpose)
      nbOReads <- nbOReads + length(alignmentInChr)
      message("Number of reads: ",length(alignmentInChr))
      reads <- keepReadCov(start(alignmentInChr),end(alignmentInChr),as.vector(strand(alignmentInChr)),keepWinPos[[chromosomeIndex]],keepWinNeg[[chromosomeIndex]],lenSeq[chromosomeIndex],win,step)
      remove(alignmentInChr)
      message("Number of kept reads: ",length(reads))
      nbKReads <- nbKReads + length(reads)
      if (length(reads)>0){
        range <- bamRange(reader,c(chromosomeIndex-1,0,end))
        bamSave(writer,range[reads,],refid=chromosomeIndex-1)
        remove(range)
      }
      remove(reads)
      if (i==length(bamfilein)){
        keepWinPos[[chromosomeIndex]] <- c(1)
        keepWinNeg[[chromosomeIndex]] <- c(1)
      }
    }
    remove(alignment)
    bamClose(reader)
    bamClose(writer)
    levelConta <- c(levelConta,(nbOReads-nbKReads)/nbOReads)
  }
  remove(keepWinPos)
  remove(keepWinNeg)
  return (levelConta)
}

