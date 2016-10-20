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
filterOneCov <- function(bamfilein,bamfileout,chromosomes=NULL,win=1000,step=100,pvalueThreshold=0.05,minR=0,maxR=0,threshold){
  #load libraries
  library(GenomicAlignments)
  library(rbamtools)
  library(Rcpp)
  library(dplyr)
  library(magrittr)
  if (!(missing(threshold))) logitThreshold <- binomial()$linkfun(threshold) 
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
    chromosomeIndex <- which(idSeq==chr)
    len <- lenSeq[chromosomeIndex]
    message("Chromosome ",chr)
    message("Length: ",len)
    #compute the normalized value of each positive/negative window to be tested by pnorm
    if (!(missing(threshold))) windows <- computeWinCov(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),len,win,step,minR,maxR,logitThreshold)
    else windows <- computeWinCovNoThreshold(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),len,win,step,minR,maxR)
    windows$Plus <- mutate(windows$Plus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% filter(pvalue<=pvalueThreshold) #compute pvalue for positive windows
    windows$Minus <- mutate(windows$Minus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% filter(pvalue<=pvalueThreshold)#compute pvalue for negative windows
    alignmentInChr <- alignment[seqnames(alignment)==chr] #get the reads in the considering chromosome
    alignment <- alignment[seqnames(alignment)!=chr] #reduce the size of alignment (for memory purpose)
    nbOReads <- nbOReads + length(alignmentInChr) 
    message("Number of reads: ",length(alignmentInChr))
    #compute the indices of reads to be kept
    keptReads <- keepReadCov(start(alignmentInChr),end(alignmentInChr),as.vector(strand(alignmentInChr)),windows$Plus$win,windows$Minus$win,len,win,step)
    remove(windows)
    remove(alignmentInChr)
    gc()
    message("Number of kept reads: ",length(keptReads))
    nbKReads <- nbKReads + length(keptReads)
    if (length(keptReads)>0){
      #get the range of kept reads
      range <- bamRange(reader,c(chromosomeIndex-1,0,len))
      #write the kept reads into output file
      bamSave(writer,range[keptReads,],refid=chromosomeIndex-1)
      remove(range)
    }
    remove(keptReads)
    gc()
  }
  bamClose(writer)
  remove(alignment)
  remove(covPos)
  remove(covNeg)
  gc()
  bamClose(reader)
  endTime2 <- proc.time()
  message("Total elapsed time ",(endTime2-startTime)[[3]]/60," minutes")
  return ((nbOReads-nbKReads)/nbOReads)
}
