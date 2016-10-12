#' filter one bamfile based on the strand specifict of itself. Strand of sliding windows are caculated based on coverage.
#' @param bamfilein the input bam file to be filterd
#' @param bamfileout the output bam file
#' @param chromosomes the list of chromosomes to be filtered
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' @param threshold the threshold upper which we keep the reads
#' @param pvalueThreshold the threshold for the p-value
#' @param minR if a window has least than minR reads then it will be all cleaned
#' @export
filterOneCov <- function(bamfilein,bamfileout,chromosomes=NULL,win=1000,step=100,threshold,pvalueThreshold=0.05,minR=0){
  #load libraries
  library(GenomicAlignments)
  library(rbamtools)
  library(Rcpp)
  library(dplyr)
  library(magrittr)
  logitThreshold <- binomial()$linkfun(threshold) 
  startTime <- proc.time()
  alignment <- readGAlignments(bamfilein) # read the input alignment
  endTime1 <- proc.time()
  message("Time loading file ",bamfilein," : ",(endTime1-startTime)[[3]]/60," minutes")
  covPos<-alignment[strand(alignment)=="+"] %>% coverage() #compute coverage of positive reads
  covNeg<-alignment[strand(alignment)=="-"] %>% coverage() #compute coverage of negative reads
  idSeq <- levels(seqnames(alignment)) #get the chromosome list
  if (is.null(chromosomes)) chromosomes <- idSeq
  lenSeq<-sapply(c(1:length(covPos)),function(i) length(covPos[[i]])) #get the length of each chromosome as the number of bases
  reader <- bamReader(bamfilein,idx=TRUE) #open a reader of the input bamfile to extract read afterward
  header <- getHeader(reader) #get the header of the input bam file
  writer <- bamWriter(header,bamfileout) #prepare to write the output bamfile with the same header
  nbOReads <- 0 #number of original reads
  nbKReads <- 0 #number of kept reads
  for (chr in chromosomes){ #filter on each chromosome
    chromosomeIndex <- which(chromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    message("Chromosome ",chr)
    message("Length: ",len)
    #compute the normalized value of each positive/negative window to be tested by pnorm
    windows <- computeWinCov(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),len,win,step,logitThreshold,minR)
    windows$Plus["pvalue"] <- pnorm(windows$Plus$value,lower.tail = FALSE) #compute pvalue for positive windows
    windows$Minus["pvalue"] <- pnorm(windows$Minus$value,lower.tail = FALSE)#compute pvalue for negative windows
    keepWinPos <- filter(windows$Plus, pvalue <= pvalueThreshold)$win # the indices of positive windows to be kept
    keepWinNeg <- filter(windows$Minus, pvalue <= pvalueThreshold)$win # the indices of negative windows to be kept
    remove(windows)
    alignmentInChr <- alignment[seqnames(alignment)==chr] #get the reads in the considering chromosome
    alignment <- alignment[seqnames(alignment)!=chr] #reduce the size of alignment (for memory purpose)
    nbOReads <- nbOReads + length(alignmentInChr) 
    message("Number of reads: ",length(alignmentInChr))
    #compute the indices of reads to be kept
    reads <- keepReadCov(start(alignmentInChr),end(alignmentInChr),as.vector(strand(alignmentInChr)),keepWinPos,keepWinNeg,lenSeq[chromosomeIndex],win,step)
    remove(keepWinPos)
    remove(keepWinNeg)
    remove(alignmentInChr)
    message("Number of kept reads: ",length(reads))
    nbKReads <- nbKReads + length(reads)
    if (length(reads)>0){
      #get the range of kept reads
      range <- bamRange(reader,c(chromosomeIndex-1,0,lenSeq[chromosomeIndex]))
      #write the kept reads into output file
      bamSave(writer,range[reads,],refid=chromosomeIndex-1)
      remove(range)
    }
    remove(reads)
  }
  bamClose(writer)
  remove(alignment)
  remove(covPos)
  remove(covNeg)
  bamClose(reader)
  endTime2 <- proc.time()
  message("Total elapsed time ",(endTime2-startTime)[[3]]/60," minutes")
  return ((nbOReads-nbKReads)/nbOReads)
}
