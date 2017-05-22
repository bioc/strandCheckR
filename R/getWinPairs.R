#' @title Get The Information of Every Window Generated From Input Bam File
#' @param bamfilein the input bam file 
#' @param histPlotFile the file name of the histogram plot
#' @param winPlotFile the file name of the windows plot
#' @param chromosomes the list of chromosomes to be considered
#' @param readLength the average length of the reads in the input bam file
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' 
#' @examples 
#' bamfilein <- system.file("data","s1.chr1.bam",package = "rnaCleanR")
#' allWin <- getWin(bamfilein,readLength = 100)
#' @export
#' 
getWinPairs <- function(bamfilein,chromosomes, readLength,win,step){
  
  alignment <- GenomicAlignments::readGAlignments(bamfilein,use.names = TRUE,param = ScanBamParam(what="flag"))
  firstReadIndex <- ((floor(alignment@elementMetadata$flag/64) %% 2) == 1) 
  alignmentFirst <- alignment[firstReadIndex]
  alignmentSecond <- alignment[!firstReadIndex]  
  
  
  covPosFirst <- alignmentFirst[strand(alignmentFirst)=="+"] %>% GenomicAlignments::coverage() 
  covNegFirst <- alignmentFirst[strand(alignmentFirst)=="-"] %>% GenomicAlignments::coverage()
  
  covPosSecond <- alignmentSecond[strand(alignmentSecond)=="+"] %>% GenomicAlignments::coverage() 
  covNegSecond <- alignmentSecond[strand(alignmentSecond)=="-"] %>% GenomicAlignments::coverage() 
  
  if (missing(win)){ 
    win <- ifelse(missing(readLength),1000,10*readLength)
  }
  if (missing(step)){
    step <- ifelse(missing(readLength),100,readLength)
  }
  if (missing(chromosomes)){
    chromosomes <- allChromosomes
  }
  
  emptyPosChromosomes1 <- sapply(covPosFirst,function(ch){length(runValue(ch))==1 && runValue(ch)[1]==0})
  emptyPosChromosomes2 <- sapply(covPosSecond,function(ch){length(runValue(ch))==1 && runValue(ch)[1]==0})
  emptyNegChromosomes1 <- sapply(covNegFirst,function(ch){length(runValue(ch))==1 && runValue(ch)[1]==0})
  emptyNegChromosomes2 <- sapply(covNegSecond,function(ch){length(runValue(ch))==1 && runValue(ch)[1]==0})
  notEmptyChromosomes <- allChromosomes[!(emptyPosChromosomes1 * emptyPosChromosomes2 * emptyNegChromosomes1 * emptyNegChromosomes2)]
  
  chromosomes <- intersect(chromosomes,notEmptyChromosomes)
  
  windowsFirst <-data.frame("Chr"=c(),"Start"=c(),"Proportion"=c(),"NbReads"=c(),"MaxCoverage"=c())
  windowsSecond <-data.frame("Chr"=c(),"Start"=c(),"Proportion"=c(),"NbReads"=c(),"MaxCoverage"=c())
  
  for (chr in chromosomes){ #walk through each chromosome
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    #compute information in each window
    windowsFirst <- rbind(windowsFirst,data.frame("Type"="First","Chr"=chr,computeWinInfo(runLength(covPosFirst[[chromosomeIndex]]),runValue(covPosFirst[[chromosomeIndex]]),runLength(covNegFirst[[chromosomeIndex]]),runValue(covNegFirst[[chromosomeIndex]]),readLength,len,win,step)))
    windowsSecond <- rbind(windowsSecond,data.frame("Type"="Second","Chr"=chr,computeWinInfo(runLength(covPosSecond[[chromosomeIndex]]),runValue(covPosSecond[[chromosomeIndex]]),runLength(covNegSecond[[chromosomeIndex]]),runValue(covNegSecond[[chromosomeIndex]]),readLength,len,win,step)))
  }
  return(rbind(windowsFirst,windowsSecond))
}
