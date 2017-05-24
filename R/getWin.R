#' @title Get The Information of Every Window Generated From an Input Single End Bam File
#' @description Return a data frame containing the information of every window generated from the input bam file. The bam file should be single end.
#' 
#' @param bamfilein the input bam file 
#' @param chromosomes the list of chromosomes to be considered
#' @param readLength the average length of the reads in the input bam file
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' 
#' @seealso getWinPairs
#' 
#' @examples 
#' bamfilein <- system.file("data","s1.chr1.bam",package = "rnaCleanR")
#' allWin <- getWin(bamfilein,readLength = 50)
#' @export
#' 
getWin <- function(bamfilein,chromosomes, readLength,win,step,minCov=0){

  # read the input alignments and compute positive/negative coverge
  alignment <- GenomicAlignments::readGAlignments(bamfilein) 
  #get the names of all chromosomes
  allChromosomes <-alignment@seqinfo@seqnames
  #get the length of each chromosome 
  lenSeq <- alignment@seqinfo@seqlengths
  
  covPos <- alignment[strand(alignment)=="+"] %>% GenomicAlignments::coverage() 
  covNeg <- alignment[strand(alignment)=="-"] %>% GenomicAlignments::coverage() 
  
  if (missing(win)){ 
    win <- ifelse(missing(readLength),1000,10*readLength)
  }
  if (missing(step)){
    step <- ifelse(missing(readLength),100,readLength)
  }
  if (missing(chromosomes)){
    chromosomes <- allChromosomes
  }
  emptyPosChromosomes <- sapply(covPos,function(ch){length(runValue(ch))==1 && runValue(ch)[1]==0})
  emptyNegChromosomes <- sapply(covNeg,function(ch){length(runValue(ch))==1 && runValue(ch)[1]==0})
  notEmptyChromosomes <- allChromosomes[!(emptyPosChromosomes * emptyNegChromosomes)]
  chromosomes <- intersect(chromosomes,notEmptyChromosomes)
  
  windows <-data.frame("Chr"=c(),"Start"=c(),"Proportion"=c(),"NbReads"=c(),"MaxCoverage"=c())
  
  
  for (chr in chromosomes){ #walk through each chromosome
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    #compute information in each window
    windows <- rbind(windows,data.frame("Chr"=chr,computeWinInfo(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),readLength,len,win,step,minCov)))
  }
  return(windows)
}
