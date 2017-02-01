#' @title Filter One Bam File 
#' 
#' @description Filter putative double strand DNA from a strand specific RNA-seq using a window sliding across the genome. 

#' 
#' @param bamfilein the input bam file to be filterd. Your bamfile should be sorted and have an index file located at the same path as well.
#' @param bamfileout the output filtered bam file
#' @param statfile the file to write the summary of the results
#' @param chromosomes the list of chromosomes to be filtered
#' @param mustKeepRanges a GRanges object defines the ranges such that every read maps to those ranges must be kept.
#' @param histPlot if TRUE then a histogram of positive proportion over all window will be generated. It's FALSE by default.
#' @param winPlot if TRUE then a plot of sum vs positive proportion over all window will be generated. It's FALSE by default.
#' @param histPlotFile the file to write the histogram plot when histPlot is TRUE
#' @param winPlotFile the file to write the window plot when winPlot is TRUE
#' @param readLength the average length of the reads
#' @param win the length of the sliding window
#' @param step the step length to sliding the window
#' @param threshold the threshold upper which we keep the reads. 0.7 by default
#' @param pvalueThreshold the threshold for the p-value. 0.05 by default
#' @param breaks the breaks of histogram when histPlot is TRUE
#' @param minCov if a window has the max coverage least than minCov, then it will be rejected regardless the strand proportion. 0 by default
#' @param maxCov if a window has the max coverage greater than maxCov, then it will be kept regardless the strand proportion. If 0 then it doesn't have effect on selecting window. 0 by default.
#' @param errorRate the probability that an RNA read takes the false strand. 0.01 by default
#' 
#' 
#' @details filterOne reads a bam file containing strand specific RNA reads, and filter putative double strand DNA. 
#' Using a window sliding across the genome, we calculate the positive/negative proportion of reads in that window.
#' Each alignment will be kept with a probability depending on the positive/negative proportions of all windows containing that alignment.
#' 
#' @seealso getPlot, filterMulti
#' 
#' @examples  
#' bamfilein <- system.file("data","s1.chr1.bam",package = "rnaCleanR")
#' filterOne(bamfilein,bamfileout="out.bam",statfile = "out.stat",histPlot = TRUE,histPlotFile = "hist.pdf", readLength = 100, threshold = 0.7)
#' @export
#' 
filterOne <- function(bamfilein,bamfileout,statfile,chromosomes,mustKeepRanges,histPlot=FALSE,winPlot=FALSE,histPlotFile,winPlotFile,readLength,win,step,pvalueThreshold=0.05,minCov=0,maxCov=0,breaks=100,threshold=0.7,errorRate=0.01){
  startTime <- proc.time()
  logitThreshold <- binomial()$linkfun(threshold) 
  if (missing(win)){ 
    win <- ifelse(missing(readLength),1000,10*readLength)
  }
  if (missing(step)){
    step <- ifelse(missing(readLength),100,readLength)
  }
  if (missing(statfile)){
    message("Summary is written to file out.stat")
    statfile <- "out.stat"
  }
  
  # read the input alignments and compute positive/negative coverge
  alignment <- GenomicAlignments::readGAlignments(bamfilein) 
  allChromosomes <-seqlevels(alignment)
  if (missing(chromosomes)){
    chromosomes <- allChromosomes
  }  
  
  covPos<-alignment[strand(alignment)=="+"] %>% GenomicAlignments::coverage() 
  covNeg<-alignment[strand(alignment)=="-"] %>% GenomicAlignments::coverage() 
  
  #get the length of each chromosome 
  lenSeq<-sapply(covPos,function(covChr) length(covChr)) 
  
  reader <- rbamtools::bamReader(bamfilein,idx=TRUE) #open a reader of the input bamfile to extract read afterward
  header <- rbamtools::getHeader(reader) #get the header of the input bam file
  writer <- rbamtools::bamWriter(header,bamfileout) #prepare to write the output bamfile with the same header
  
  nbOReads <- 0 #number of original reads
  nbKReads <- 0 #number of kept reads
  
  append=FALSE 
  allWin <- data.frame("chr"=c(),"propor"=c(),"sum"=c(),"max"=c(),"group"=c()) #data frame contains information of all windows for plotting
  mustKeepPosWin <- c()
  mustKeepNegWin <- c()
  for (chr in chromosomes){ #walk through each chromosome
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    #compute strand information in each window
    if (winPlot==TRUE || histPlot==TRUE){#get details of each window for the plots
      windows <- computeWinVerbose(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),readLength,len,win,step,minCov,maxCov,logitThreshold)
      windows$Plus <- dplyr::mutate(windows$Plus,"chr"=chr)
      windows$Minus <- dplyr::mutate(windows$Minus,"chr"=chr)
      dupWin <- duplicated(c(windows$Plus$win,windows$Minus$win))
      allWin <- rbind(allWin,rbind(dplyr::select(windows$Plus,c(chr,propor,sum,max,group)),dplyr::select(windows$Minus,c(chr,propor,sum,max,group)))[!dupWin,])
    }
    else{
      windows <- computeWin(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),readLength,len,win,step,minCov,maxCov,logitThreshold)
    }
    if (!missing(mustKeepRanges)){
      mkr <- mustKeepRanges[seqnames(mustKeepRanges)==chr]
      rg <- ranges(mkr)
      pos <- which(strand(mkr)=="+")
      mustKeepPosWin <- getWin(rg[pos,],win,step)  
      mustKeepNegWin <- getWin(rg[-pos,],win,step)  
      remove(mkr)
      remove(rg)
      remove(pos)
    }
    windows$Plus <- dplyr::mutate(windows$Plus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% dplyr::filter(pvalue<=pvalueThreshold | win %in% mustKeepPosWin) %>% dplyr::select(c(win,propor))
    windows$Minus <- dplyr::mutate(windows$Minus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% dplyr::filter(pvalue<=pvalueThreshold | win %in% mustKeepNegWin) %>% dplyr::select(c(win,propor))
    
    keptReads <- c()
    nbOReadsChr <- 0
    if (nrow(windows$Plus)>0 || nrow(windows$Minus)>0){
      alignmentInChr <- alignment[seqnames(alignment)==chr] #get the reads in the considering chromosome
      alignment <- alignment[seqnames(alignment)!=chr] #reduce the size of alignment, for memory purpose
      nbOReadsChr <- length(alignmentInChr)
      nbOReads <- nbOReads + nbOReadsChr
      
      strand <- as.vector(strand(alignmentInChr))
      indexPos <- which(strand=="+")
      indexNeg <- which(strand=="-")
      remove(strand)
      fragments <- getFragment(alignmentInChr)
      remove(alignmentInChr)
      
      keptReads <- keepRead(fragments$Pos,fragments$Neg,windows$Plus,windows$Minus,win,step,errorRate);   
      remove(windows)
      keptReads <- c(indexPos[unique(keptReads$Pos)],indexNeg[unique(keptReads$Neg)]) %>% sort() 
      cat(paste0("Sequence ",chr,", length: ",len,", number of reads: ",nbOReadsChr,", number of kept reads: ",length(keptReads),"\n"),file=statfile,append=append)
      if (append==FALSE) append <- TRUE  
    }
    if (length(keptReads)>0){##write the kept reads into output file
      #get the range of kept reads
      nbKReads <- nbKReads + length(keptReads)
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
