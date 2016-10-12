#' filter several bamfiles based on the strand specific of the merged files. Strand of sliding windows are caculated based on coverage.
#' @param bamfilein the input bam files to be filterd
#' @param bamfileout the output bam files
#' @param chromosomes the list of chromosomes to be filtered
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' @param threshold the threshold upper which we keep the reads
#' @param pvalueThreshold the threshold for the p-value
#' @export
filterMultiCov <- function(bamfilein,bamfileout,chromosomes=NULL,win=1000,step=100,threshold,pvalueThreshold=0.05,minR=0){
  library(GenomicAlignments)
  library(rbamtools)
  library(Rcpp)
  library(dplyr)
  #sourceCpp("computeWinCov.cpp")
  #sourceCpp("keepRead.cpp")
  logitThreshold <- binomial()$linkfun(threshold)
  alignment <- readGAlignments(bamfilein[1])
  covPos <- alignment[strand(alignment)=="+"] %>% coverage()
  covNeg <- alignment[strand(alignment)=="-"] %>% coverage()
  idSeq <- levels(seqnames(alignment))
  remove(alignment)
  if (is.null(chromosomes)) chromosomes <- idSeq
  lenSeq<-sapply(c(1:length(covPos)),function(i) length(covPos[[i]]))
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
    chromosomeIndex <- which(chromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    windows <- computeWinCov(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),len,win,step,logitThreshold,minR)
    windows$Plus["pvalue"] <- pnorm(windows$Plus$value,lower.tail = FALSE)
    windows$Minus["pvalue"] <- pnorm(windows$Minus$value,lower.tail = FALSE)
    keepWinPos[[chromosomeIndex]] <- filter(windows$Plus, pvalue <= pvalueThreshold)$win
    keepWinNeg[[chromosomeIndex]] <- filter(windows$Minus, pvalue <= pvalueThreshold)$win
    remove(windows)
  }
  remove(covPos)
  remove(covNeg)
  levelConta <- c()
  #write the kept reads to output files
  for (i in c(1:length(bamfilein))){
    message("File ",bamfilein[i])
    alignment <- readGAlignments(bamfilein[i])
    reader <- bamReader(bamfilein[i],idx=TRUE)
    header <- getHeader(reader)
    writer <- bamWriter(header,bamfileout[i])
    nbOReads <- 0
    nbKReads <- 0
    for (chr in chromosomes){
      message("Chromosome ",chr)
      chromosomeIndex <- which(chromosomes==chr)
      end <- lenSeq[chromosomeIndex]
      alignmentInChr <- alignment[seqnames(alignment)==chr]
      alignment <- alignment[seqnames(alignment)!=chr]
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

