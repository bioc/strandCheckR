keepCount <- function(bamfilein,bamfileout,chromosomes,allChromosomes,lenSeq,win,step,threshold,pvalueThreshold,minR,maxR,limit){#compute the reads to be kept for the whole genome
  library(GenomicAlignments)
  library(Rcpp)
  library(dplyr)
  library(magrittr)
  if (!missing(threshold)) logitThreshold <- binomial()$linkfun(threshold)
  alignments <- list()
  keepReads <- list() # the list of kept reads
  for (i in c(1:length(bamfilein))){#load bam file
    alignments[[i]] <- readGAlignments(bamfilein[i],param=ScanBamParam(what=c("cigar")))
    keepReads[[i]] <- list()
  }
  for (chr in chromosomes){
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    kch <- keepCountChr(length(bamfilein),alignments,chr,len,win,step,pvalueThreshold,minR,maxR,limit,logitThreshold)
    for (i in c(1:length(bamfilein))){
      keepReads[[i]][[chr]] <- kch[[i]]
      kch[[i]] <-c(1)
    }
    remove(kch)
  }
  rm(list=setdiff(ls(), "keepReads"))
  gc()
  return (keepReads)
}