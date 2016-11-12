#' @title Calculate the Reads to be Kept
#' 
#' @description None
#' 
#' @param bamfilein the input bam files to be filterd
#' @param chromosomes the list of chromosomes to be filtered
#' @param allChromosomes the list of all chromosomes
#' @param lenSeq the length of every chromosome
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' @param pvalueThreshold the threshold for the p-value
#' @param minCov if a window has the max coverage least than minCov, then it will be rejected
#' @param maxCov if a window has the max coverage greater than maxCov, then it will be kept
#' @param limit the proportion that a read must be inside a window in order to be counted
#' @param threshold the threshold upper which we keep the reads
#' 

keepCount <- function(bamfilein,chromosomes,allChromosomes,lenSeq,win,step,pvalueThreshold,minCov,maxCov,limit,threshold){#compute the reads to be kept for the whole genome
  keepReads <- list() # the list of kept reads
  mergedAlignments <- list() # the list of input reads
  for (chr in chromosomes){
    mergedAlignments[[chr]] <- list() 
  }
  for (i in seq_along(bamfilein)){#load bam file
    al <- readGAlignments(bamfilein[i],param=ScanBamParam(what=c("cigar")))
    for (chr in chromosomes){
      mergedAlignments[[chr]][[i]] <- al[seqnames(al)==chr]
    }
    remove(al)
    keepReads[[i]] <- list()
  }
  for (chr in chromosomes){
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    message("Chromosome: ",chr,", Length: ",len)
    keepReadChr <- keepCountChr(length(bamfilein),mergedAlignments[[chr]],len,win,step,pvalueThreshold,minCov,maxCov,limit,threshold)
    mergedAlignments[[chr]] <- c()
    for (i in c(1:length(bamfilein))){
      keepReads[[i]][[chr]] <- keepReadChr[[i]]
      keepReadChr[[i]] <-c(1)
    }
    remove(keepReadChr)
  }
  rm(list=setdiff(ls(), "keepReads"))
  gc()
  return (keepReads)
}