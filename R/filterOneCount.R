#' filter one bamfile based on the strand specifict of itself. Strand of sliding windows are caculated based on read counts.
#' @param bamfilein the input bam file to be filterd
#' @param bamfileout the output bam file
#' @param chromosomes the list of chromosomes to be filtered
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' @param threshold the threshold upper which we keep the reads
#' @param pvalueThreshold the threshold for the p-value
#' @param minR if a window has least than minR reads then it will be all cleaned
#' @param limit the proportion of a read that it should not exceed to be considered to be in a window
#' @export
filterOneCount <- function(bamfilein,bamfileout,chromosomes=NULL,win=1000,step=100,threshold,pvalueThreshold=0.05,minR=0,maxR=0,limit=0.25){
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
  writer <- bamWriter(header,bamfileout) #prepare to write the output bamfile with the same headermessage("Chromosome ",chr)
  for (chr in chromosomes){
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    startTime <- proc.time()
    position <- computePosition(alignment,chr)
    if (minR==0 && maxR==0) {
      windows <- computeWinCount0(position$Pos$start,position$Pos$end,position$Neg$start,position$Neg$end,len,win,step,logitThreshold,limit)
    }
    else {
      position <- reorder(position,win,step,limit)
      windows <- computeWinCount(position$Pos$start,position$Pos$end,position$Neg$start,position$Neg$end,len,win,step,logitThreshold,minR,maxR,limit)
    }
    windows$Plus <- mutate(windows$Plus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% filter(pvalue<=pvalueThreshold) #compute pvalue for positive windows
    windows$Minus <- mutate(windows$Minus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% filter(pvalue<= pvalueThreshold) #compute pvalue for negative windows
    #message("keep win Pos: ",nrow(windows$Plus))
    #message("keep win Neg: ",nrow(windows$Minus))
    #compute the indices of reads to be kept
    keptFrags <- keepReadOne(position$Pos$start,position$Pos$end,as.integer(position$Pos$group),position$Neg$start,position$Neg$end,as.integer(position$Neg$group),windows$Plus$win,windows$Minus$win,len,win,step,limit)
    remove(windows)
    gc()
    #message("keep frag Pos: ",length(keptFrags$Pos))
    #message("keep frag Neg: ",length(keptFrags$Neg))
    keptReads <- c(position$Index$Pos[unique(keptFrags$Pos)],position$Index$Neg[unique(keptFrags$Neg)]) %>% sort()
    remove(keptFrags)
    gc()
    message("Chromosome: ",chr, ", Length: ",len, ", Number of reads: ",position$nbRead,", Number of kept reads: ",length(keptReads))
    remove(position)
    gc()
    if (length(keptReads)>0){
      #get the range of kept reads
      range <- bamRange(reader,c(chromosomeIndex-1,0,len))
      #write the kept reads into output file
      bamSave(writer,range[keptReads,],refid=chromosomeIndex-1)
      remove(range)
    }
    remove(keptReads)
  }
  remove(alignment)
  bamClose(writer)
  bamClose(reader)
  endTime2 <- proc.time()
  message("Total elapsed time ",(endTime2-startTime)[[3]]/60," minutes")
}
