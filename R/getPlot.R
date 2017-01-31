#' @title Get The Plots Generated From Input Bam File
#' 
#' @description Generate the histogram of the positive proportions, and the plot presenting number of reads vs positive proportion over all sliding windows.

#' @param bamfilein the input bam file 
#' @param histPlotFile the file name of the histogram plot
#' @param winPlotFile the file name of the windows plot
#' @param chromosomes the list of chromosomes to be considered
#' @param readLength the average length of the reads in the input bam file
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' @param breaks the breaks of the histogram plot
#' @param xlim xlim of the window plot
#' @param minCov the coverage of a window under which the window is ignored
#' 
#' @examples 
#' bamfilein <- system.file("data","s1.chr1.bam",package = "rnaCleanR")
#' getPlot(bamfilein,histPlotFile = "hist.pdf", winPlotFile = "win.pdf",readLength = 100)
#' @export
#' 
getPlot <- function(bamfilein,histPlotFile,winPlotFile,chromosomes = NULL, readLength,win,step,breaks = 100,xlim,minCov=0){
  # read the input alignments and compute positive/negative coverge
  if (missing(win)){ 
    win <- ifelse(missing(readLength),1000,10*readLength)
  }
  if (missing(step)){
    step <- ifelse(missing(readLength),100,readLength)
  }
  
  
  alignment <- GenomicAlignments::readGAlignments(bamfilein) 
  covPos <- alignment[strand(alignment)=="+"] %>% GenomicAlignments::coverage() 
  covNeg <- alignment[strand(alignment)=="-"] %>% GenomicAlignments::coverage() 
  
  #get the names of all chromosomes
  allChromosomes <- levels(seqnames(alignment)) 
  if (is.null(chromosomes)) chromosomes <- allChromosomes
  
  #get the length of each chromosome 
  lenSeq<-sapply(covPos,function(covChr) length(covChr)) 
  
  windows <-data.frame("win"=c(), "chr"=c(),"propor"=c(),"sum"=c(),"max"=c(),"group"=c())
  
  for (chr in chromosomes){ #filter on each chromosome
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    #compute strand information in each window
    windows <- rbind(windows,computeWinPlot(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),readLength,len,win,step,minCov) %>% dplyr::mutate("chr"=chr))
  }
  if (missing(histPlotFile)){
    histPlotFile <- paste0(bamfilein,"_hist.pdf")
  }
  histPlot(windows,histPlotFile,breaks = breaks)
  if (missing(winPlotFile)){
    winPlotFile <- paste0(bamfilein,"_win.pdf")
  }
  if (missing(xlim)){
    winPlot(windows,winPlotFile)
  }
  else {
    winPlot(windows,winPlotFile,xlim)
  }
  remove(alignment)
  remove(covPos)
  remove(covNeg)
}
