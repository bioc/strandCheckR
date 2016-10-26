keepCount <- function(bamfilein,bamfileout,chromosomes,allChromosomes,lenSeq,win,step,pvalueThreshold,minR,maxR,limit,threshold){#compute the reads to be kept for the whole genome
  library(GenomicAlignments)
  library(Rcpp)
  library(dplyr)
  library(magrittr)
  alignments <- list()
  keepReads <- list() # the list of kept reads
  alignments <- list() # the list of input reads
  for (chr in chromosomes){
    alignments[[chr]] <- list() 
  }
  for (i in c(1:length(bamfilein))){#load bam file
    al <- readGAlignments(bamfilein[i],param=ScanBamParam(what=c("cigar")))
    for (chr in chromosomes){
      alignments[[chr]][[i]] <- al[seqnames(al)==chr]
    }
    remove(al)
    keepReads[[i]] <- list()
  }
  for (chr in chromosomes){
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    message("Chromosome: ",chr,", Length: ",len)
    kch <- keepCountChr(length(bamfilein),alignments[[chr]],len,win,step,pvalueThreshold,minR,maxR,limit,threshold)
    alignments[[chr]] <- c()
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