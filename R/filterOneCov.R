#' @title Filter One Bam File Based On Coverage
#' 
#' @description filter one bamfile based on the strand specifict of itself. 
#' Strand of sliding windows are caculated based on the sum of coverage.
#' 
#' @details None
#' 
#' @param bamfilein the input bam file to be filterd
#' @param bamfileout the output bam file
#' @param chromosomes the list of chromosomes to be filtered
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' @param threshold the threshold upper which we keep the reads
#' @param pvalueThreshold the threshold for the p-value
#' @param minCov if a window has the max coverage least than minCov, then it will be rejected
#' @param maxCov if a window has the max coverage greater than maxCov, then it will be kept
#' 
#' @importFrom GenomicAlignments readGAlignments coverage 
#' @import rbamtools
#' @import Rcpp
#' @import dplyr
#' @import magrittr
#' 
#' @export
filterOneCov <- function(bamfilein,bamfileout,chromosomes=NULL,win=1000,step=100,pvalueThreshold=0.05,minCov=0,maxCov=0,threshold=0.7){
  startTime <- proc.time()
  
  # read the input alignments and compute positive/negative coverge
  alignment <- readGAlignments(bamfilein) 
  covPos<-alignment[strand(alignment)=="+"] %>% coverage() 
  covNeg<-alignment[strand(alignment)=="-"] %>% coverage() 
  
  #get the names of all chromosomes
  allChromosomes <- levels(seqnames(alignment)) 
  if (is.null(chromosomes)) chromosomes <- allChromosomes
  
  #get the length of each chromosome 
  lenSeq<-sapply(covPos,function(covChr) length(covChr)) 
  
  reader <- bamReader(bamfilein,idx=TRUE) #open a reader of the input bamfile to extract read afterward
  header <- getHeader(reader) #get the header of the input bam file
  writer <- bamWriter(header,bamfileout) #prepare to write the output bamfile with the same header
  
  nbOReads <- 0 #number of original reads
  nbKReads <- 0 #number of kept reads
  logitThreshold <- binomial()$linkfun(threshold) 
  
  for (chr in chromosomes){ #filter on each chromosome
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    message("Chromosome ",chr,", Length: ",len)
    
    #compute strand information in each window
    windows <- computeWinCov(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),len,win,step,minCov,maxCov,logitThreshold)
    windows$Plus <- mutate(windows$Plus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% filter(pvalue<=pvalueThreshold) 
    windows$Minus <- mutate(windows$Minus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% filter(pvalue<=pvalueThreshold)
    
    alignmentInChr <- alignment[seqnames(alignment)==chr] #get the reads in the considering chromosome
    alignment <- alignment[seqnames(alignment)!=chr] #reduce the size of alignment (for memory purpose)
    nbOReads <- nbOReads + length(alignmentInChr) 
    message("Number of reads: ",length(alignmentInChr))
    
    keptReads <- keepReadCov(start(alignmentInChr),end(alignmentInChr),as.vector(strand(alignmentInChr)),windows$Plus$win,windows$Minus$win,len,win,step)#compute the indices of reads to be kept
    remove(windows)
    remove(alignmentInChr)
    message("Number of kept reads: ",length(keptReads))
    nbKReads <- nbKReads + length(keptReads)
    
    if (length(keptReads)>0){##write the kept reads into output file
      #get the range of kept reads
      range <- bamRange(reader,c(chromosomeIndex-1,0,len))
      #write the kept reads into output file
      bamSave(writer,range[keptReads,],refid=chromosomeIndex-1)
      remove(range)
    }
    remove(keptReads)
  }
  bamClose(writer)
  remove(alignment)
  remove(covPos)
  remove(covNeg)
  bamClose(reader)
  endTime2 <- proc.time()
  message("Total elapsed time ",(endTime2-startTime)[[3]]/60," minutes")
  return ((nbOReads-nbKReads)/nbOReads) #return proportion of removed reads
}
