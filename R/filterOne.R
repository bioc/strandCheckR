#' @title Filter One Bam File Based On Coverage
#' 
#' @description Filter double strand DNA from a strand specific RNA-seq. Positive/negative strands of sliding windows are caculated based on the coverage of positive/negative reads that overlap the window.

#' 
#' @param bamfilein the input bam file to be filterd. Your bamfile should be sorted and have an index file located at the same path as well.
#' @param bamfileout the output bam file
#' @param statfile the file to write the summary of the results
#' @param histfile the file to write the histogram of the positive proportions over all windows
#' @param chromosomes the list of chromosomes to be filtered
#' @param win the length of the sliding window. 1000 by default
#' @param step the step length to slide the window. 100 by default
#' @param threshold the threshold upper which we keep the reads. 0.7 by default
#' @param pvalueThreshold the threshold for the p-value. 0.05 by default
#' @param breaks the breaks to plot histogram
#' @param minCov if a window has the max coverage least than minCov, then it will be rejected regardless the strand proportion. 0 by default
#' @param maxCov if a window has the max coverage greater than maxCov, then it will be kept regardless the strand proportion. If 0 then it doesn't have effect on selecting window. 0 by default.
#' 
#' 
#' @details filterOneCov reads a bam file containing strand specific RNA reads, and filter the potential double strand contamination DNA from that. 
#' Using a window sliding across the genome, we estimate the probability to keep or reject it based on the proportion of positive/negative reads in the window.
#' 
#' @import Rcpp
#' @import magrittr
#' @import S4Vectors
#' @import BiocGenerics
#' 
#' @examples  
#' bamfilein <- system.file("data","chr1.bam",package = "rnaCleanR")
#' filterOne(bamfilein,bamfileout="out.bam",statfile = "out.stat",histPlot = TRUE,histfile = "hist.png", readLength = 100, threshold = 0.7)
#' @export
#' 
filterOne <- function(bamfilein,bamfileout,statfile,histPlot=FALSE,winPlot=FALSE,histPlotFile,winPlotFile,chromosomes=NULL,readLength,win,step,pvalueThreshold=0.05,minCov=0,maxCov=0,breaks=100,threshold=0.7,errorRate=0.01){
  startTime <- proc.time()
  if (missing(win)){ 
    if (missing(readLength)){ 
      win = 1000
    }
    else {
      win = 10*readLength
    }
  }
  if (missing(step)){
    if (missing(readLength)){
      step = 100
    }
    else{
      step = readLength  
    }
  }
  if (missing(statfile)){
    message("Summary is written to file out.stat")
    statfile <- "out.stat"
  }
  
  # read the input alignments and compute positive/negative coverge
  alignment <- GenomicAlignments::readGAlignments(bamfilein) 
  covPos<-alignment[strand(alignment)=="+"] %>% GenomicAlignments::coverage() 
  covNeg<-alignment[strand(alignment)=="-"] %>% GenomicAlignments::coverage() 
  
  #get the names of all chromosomes
  allChromosomes <- levels(seqnames(alignment)) 
  if (is.null(chromosomes)) chromosomes <- allChromosomes
  
  #get the length of each chromosome 
  lenSeq<-sapply(covPos,function(covChr) length(covChr)) 
  
  reader <- rbamtools::bamReader(bamfilein,idx=TRUE) #open a reader of the input bamfile to extract read afterward
  header <- rbamtools::getHeader(reader) #get the header of the input bam file
  writer <- rbamtools::bamWriter(header,bamfileout) #prepare to write the output bamfile with the same header
  
  nbOReads <- 0 #number of original reads
  nbKReads <- 0 #number of kept reads
  logitThreshold <- binomial()$linkfun(threshold) 
  append=FALSE
  allWin <- data.frame("win"=c(),"propor"=c(),"sum"=c(),"max"=c(),"group"=c())
  for (chr in chromosomes){ #filter on each chromosome
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    #compute strand information in each window
    if (histPlot==TRUE || winPlot==TRUE){#get details of each window for the plots
      windows <- computeWin1(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),readLength,len,win,step,minCov,maxCov,logitThreshold)
      dupWin <- duplicated(c(windows$Plus$win,windows$Minus$win))
      allWin <- rbind(allWin,rbind(dplyr::select(windows$Plus,c(win,propor,sum,max,group)),dplyr::select(windows$Minus,c(win,propor,sum,max,group)))[!dupWin,])
    }
    else{
      windows <- computeWin(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),readLength,len,win,step,minCov,maxCov,logitThreshold)
    }
    maxWin <- 0
    if (nrow(windows$Plus)>0) maxWin <- windows$Plus$win[nrow(windows$Plus)]
    if (nrow(windows$Minus)>0) maxWin <- max(maxWin,windows$Minus$win[nrow(windows$Minus)])
    
    windows$Plus <- dplyr::mutate(windows$Plus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% dplyr::filter(pvalue<=pvalueThreshold) %>% dplyr::select(c(win,propor))
    windows$Minus <- dplyr::mutate(windows$Minus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% dplyr::filter(pvalue<=pvalueThreshold) %>% dplyr::select(c(win,propor))
    
    alignmentInChr <- alignment[seqnames(alignment)==chr] #get the reads in the considering chromosome
    alignment <- alignment[seqnames(alignment)!=chr] #reduce the size of alignment (for memory purpose)
    nbOReadsChr <- length(alignmentInChr)
    nbOReads <- nbOReads + nbOReadsChr
    
    position <- data.frame("group"=c(),"start"=c(),"end"=c(),"sample"=c()) #position of each fragment 
    alPos <- alignmentInChr[strand(alignmentInChr)=='+']
    alNeg <- alignmentInChr[strand(alignmentInChr)=='-']
    index <- getIndex(as.vector(strand(alignmentInChr)))
    remove(alignmentInChr)
    positionPos <- GenomicAlignments::extractAlignmentRangesOnReference(GenomicAlignments::cigar(alPos),pos=start(alPos)) %>% data.frame() %>% dplyr::select(-c(group_name,width))  
    positionPos <- positionPos[order(positionPos$start),]
    positionNeg <- GenomicAlignments::extractAlignmentRangesOnReference(GenomicAlignments::cigar(alNeg),pos=start(alNeg)) %>% data.frame() %>% dplyr::select(-c(group_name,width))  
    positionNeg <- positionNeg[order(positionNeg$start),]
    
    keptFrags <- c()
    if (maxWin>0){
      keptFrags <- keepRead(positionPos$start,positionPos$end,positionPos$group,positionNeg$start,positionNeg$end,positionNeg$group,windows$Plus$propor,windows$Minus$propor,windows$Plus$win,windows$Minus$win,maxWin,win,step,errorRate);   
    }
    keptReads <- c(index$Pos[unique(keptFrags$Pos)],index$Neg[unique(keptFrags$Neg)]) %>% sort() 
    remove(windows)
    if (chromosomeIndex!=1){
      append=TRUE
    }
    cat(paste0("Chromosome ",chr,", length: ",len,", number of reads: ",nbOReadsChr,", number of kept reads: ",length(keptReads),"\n"),file=statfile,append=append)
    nbKReads <- nbKReads + length(keptReads)
    
    if (length(keptReads)>0){##write the kept reads into output file
      #get the range of kept reads
      range <- rbamtools::bamRange(reader,c(chromosomeIndex-1,0,len))
      #write the kept reads into output file
      rbamtools::bamSave(writer,range[keptReads,],refid=chromosomeIndex-1)
      remove(range)
    }
    remove(keptReads)
  }
  rbamtools::bamClose(writer)
  remove(alignment)
  remove(covPos)
  remove(covNeg)
  rbamtools::bamClose(reader)
  
  cat("Summary:\n",file = statfile, append = append)
  cat(paste0("Number of original reads: ",nbOReads,", number of kept reads: ",nbKReads,", removal proportion: ",(nbOReads-nbKReads)/nbOReads),"\n",file = statfile,append=TRUE)
  if (histPlot == TRUE){
    if (missing(histPlotFile)){
      histPlotFile <- paste0(bamfilein,"_hist.pdf")
    }
    histPlot(allWin,histPlotFile,breaks = breaks)
  }
  if (winPlot == TRUE){
    if (missing(winPlotFile)){
      winPlotFile <- paste0(bamfilein,"_win.pdf")
    }
    winPlot(allWin,winPlotFile)
  }
  endTime2 <- proc.time()
  cat(paste0("Total elapsed time ",(endTime2-startTime)[[3]]/60," minutes\n"),file = statfile,append=TRUE)
}
