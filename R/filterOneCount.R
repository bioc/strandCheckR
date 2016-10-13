#' filter one bamfile based on the strand specifict of itself. Strand of sliding windows are caculated based on read counts.
#' @param bamfilein the input bam file to be filterd
#' @param bamfileout the output bam file
#' @param chromosomes the list of chromosomes to be filtered
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' @param threshold the threshold upper which we keep the reads
#' @param pvalueThreshold the threshold for the p-value
#' @param minR if a window has least than minR reads then it will be all cleaned
#' @export
filterOneCount <- function(bamfilein,bamfileout,chromosomes=NULL,win=1000,step=100,threshold,pvalueThreshold=0.05,minR=0,limit=0.25){
  #load libraries
  library(GenomicAlignments)
  library(rbamtools)
  library(Rcpp)
  library(dplyr)
  library(magrittr)
  logitThreshold <- binomial()$linkfun(threshold) 
  alignment <- readGAlignments(bamfilein,param=ScanBamParam(what=c("cigar")))
  reader <- bamReader(bamfilein,idx=TRUE) #open a reader of the input bamfile to extract read afterward
  header <- getHeader(reader) #get the header of the input bam file
  refSeqs <- getRefData(reader)
  allChromosomes <- refSeqs$SN
  if (is.null(chromosomes)) chromosomes<-allChromosomes
  lenSeq <- refSeqs$LN
  writer <- bamWriter(header,bamfileout) #prepare to write the output bamfile with the same header
  nbOReads <- 0 #number of original reads
  nbKReads <- 0 #number of kept reads
  for (chr in chromosomes){
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    message("Chromosome ",chr)
    message("Length: ",len)
    startTime <- proc.time()
    alignmentInChr <- alignment[seqnames(alignment)==chr] #get the reads in the considering chromosome
    index <- getIndex(as.vector(strand(alignmentInChr))) #get index of positive and negative reads
    nbOReads <- nbOReads + length(alignmentInChr) 
    message("Number of reads: ",length(alignmentInChr))
    #alignment <- alignment[seqnames(alignment)!=chr] #reduce the size of alignment (for memory purpose)
    alignmentInChrPos<-alignmentInChr[strand(alignmentInChr)=="+"]
    positionPos<-extractAlignmentRangesOnReference(cigar(alignmentInChrPos),pos=start(alignmentInChrPos)) %>% data.frame() %>% select(-c(group_name,width))
    remove(alignmentInChrPos)
    alignmentInChrNeg<-alignmentInChr[strand(alignmentInChr)=="-"]
    positionNeg<-extractAlignmentRangesOnReference(cigar(alignmentInChrNeg),pos=start(alignmentInChrNeg)) %>% data.frame() %>% select(-c(group_name,width))
    remove(alignmentInChrNeg)
    if (minR==0) {
      windows <- computeWinCount0(positionPos$start,positionPos$end,positionNeg$start,positionNeg$end,len,win,step,logitThreshold,limit)
    }
    else {
      positionPos[["firstW"]] <- ceiling((positionPos$start-win-1+floor((positionPos$end-positionPos$start+1)*limit))/step)
      positionNeg[["firstW"]] <- ceiling((positionNeg$start-win-1+floor((positionNeg$end-positionNeg$start+1)*limit))/step)
      positionPos <-  positionPos[order(positionPos$firstW),] #reorder positionPos following the starting position of each fragment
      positionNeg <-  positionNeg[order(positionNeg$firstW),] #reorder positionNeg following the starting position of each fragment
      windows <- computeWinCount(positionPos$start,positionPos$end,positionNeg$start,positionNeg$end,len,win,step,logitThreshold,minR,limit)
    }
    windows$Plus["pvalue"] <- pnorm(windows$Plus$value,lower.tail = FALSE) #compute pvalue for positive windows
    windows$Minus["pvalue"] <- pnorm(windows$Minus$value,lower.tail = FALSE)#compute pvalue for negative windows
    keepWinPos <- filter(windows$Plus, pvalue <= pvalueThreshold)$win # the indices of positive windows to be kept
    keepWinNeg <- filter(windows$Minus, pvalue <= pvalueThreshold)$win # the indices of negative windows to be kept
    remove(windows)
    #compute the indices of reads to be kept
    reads <- keepReadOne(positionPos$start,positionPos$end,as.integer(positionPos$group),positionNeg$start,positionNeg$end,as.integer(positionNeg$group),keepWinPos,keepWinNeg,lenSeq[chromosomeIndex],win,step,limit)
    remove(positionPos)
    remove(positionNeg)
    remove(keepWinPos)
    remove(keepWinNeg)
    reads$Pos <- unique(reads$Pos)
    reads$Neg <- unique(reads$Neg)
    reads <- c(index$Pos[reads$Pos],index$Neg[reads$Neg]) %>% sort()
    remove(index)
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
  bamClose(reader)
  endTime2 <- proc.time()
  message("Total elapsed time ",(endTime2-startTime)[[3]]/60," minutes")
  return ((nbOReads-nbKReads)/nbOReads)
}
