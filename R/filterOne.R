#' @title Filter One Bam File 
#' 
#' @description Filter putative double strand DNA from a strand specific RNA-seq using a window sliding across the genome. 

#' 
#' @param bamfilein the input bam file to be filterd. Your bamfile should be sorted and have an index file located at the same path as well.
#' @param bamfileout the output filtered bam file
#' @param statfile the file to write the summary of the results
#' @param chromosomes the list of chromosomes to be filtered
#' @param mustKeepRanges a GRanges object defines the ranges such that every read maps to those ranges must be always kept regardless the strand proportion of the windows containing them.
#' @param getWin if TRUE, the function will return a data frame containing the information of all windows as will. It's FALSE by default.
#' @param readLength the average length of the reads
#' @param win the length of the sliding window
#' @param step the step length to sliding the window
#' @param threshold the threshold upper which we keep the reads. 0.7 by default
#' @param pvalueThreshold the threshold for the p-value. 0.05 by default
#' @param minCov if a window has the max coverage least than minCov, then it will be rejected regardless the strand proportion. 0 by default
#' @param maxCov if a window has the max coverage greater than maxCov, then it will be kept regardless the strand proportion. If 0 then it doesn't have effect on selecting window. 0 by default.
#' @param errorRate the probability that an RNA read takes the false strand. 0.01 by default
#' 
#' 
#' @details filterOne reads a bam file containing strand specific RNA reads, and filter putative double strand DNA. 
#' Using a window sliding across the genome, we calculate the positive/negative proportion of reads in that window.
#' For each window, we use logistic regression to estimate the proportion of reads in the window derived from 
#' stranded RNA (positive or negative). 
#' 
#' \eqn{\pi}: proportion of reads in the window derived from stranded RNA (positive or negative)
#' 
#' Null hypothesis: \eqn{\pi < {\pi}_{0}} where \eqn{{\pi}_{0}} is the given threshold.
#' 
#' Only windows with p-value <= 0.05 are kept. Considering a positive window that is kept, let P be its number of positive reads, and let M
#' be its number of negative reads. Since these M negative reads should come from double-strand DNA, then there should be also M postive reads among the
#' P positive reads come from double-strand DNA. In other words, there are only (P-M) positive reads come from RNA. Hence, each positive read in this window is kept 
#' with the probability equalling (P-M)/P. Each negative read is kept with the probability equalling the given errorRate which is the rate that
#' an RNA read of your sample has wrong strand.
#' 
#' Since each alignment can be belonged to several windows, then the probability of keeping an alignment is the maximum probability defined by
#' all windows that contain it.
#' 
#' @seealso getPlot, filterMulti
#' 
#' @examples  
#' bamfilein <- system.file("data","s1.chr1.bam",package = "rnaCleanR")
#' filterOne(bamfilein,bamfileout="out.bam",statfile = "out.stat",readLength = 100, threshold = 0.7)
#' @export
#' 
filterOne <- function(bamfilein,bamfileout,statfile,chromosomes,mustKeepRanges,getWin=FALSE,readLength,win,step,pvalueThreshold=0.05,minCov=0,maxCov=0,threshold=0.7,errorRate=0.01){
  startTime <- proc.time() 
  
  # read the info of the input alignments and compute positive/negative coverge
  alignment <- GenomicAlignments::readGAlignments(bamfilein,use.names = TRUE) 
  allChromosomes <-alignment@seqinfo@seqnames #get the name of each chromosome 
  lenSeq <- alignment@seqinfo@seqlengths #get the length of each chromosome 
  covPos<-alignment[strand(alignment)=="+"] %>% GenomicAlignments::coverage() #calculate coverage came from positive reads 
  covNeg<-alignment[strand(alignment)=="-"] %>% GenomicAlignments::coverage() #calculate coverage came from negative reads
  
  # treat missing parameters
  if (missing(chromosomes)){
    chromosomes <- allChromosomes
  }  
  emptyPosChromosomes <- sapply(covPos,function(ch){length(runValue(ch))==1 && runValue(ch)[1]==0})
  emptyNegChromosomes <- sapply(covNeg,function(ch){length(runValue(ch))==1 && runValue(ch)[1]==0})
  notEmptyChromosomes <- allChromosomes[!(emptyPosChromosomes * emptyNegChromosomes)]
  chromosomes <- intersect(chromosomes,notEmptyChromosomes)
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
  
  reader <- rbamtools::bamReader(bamfilein,idx=TRUE) #open a reader of the input bamfile to extract read afterward
  header <- rbamtools::getHeader(reader) #get the header of the input bam file
  writer <- rbamtools::bamWriter(header,bamfileout) #prepare to write the output bamfile with the same header
  logitThresholdP <- binomial()$linkfun(threshold) 
  logitThresholdM <- binomial()$linkfun(1-threshold) 
  nbOReads <- 0 #number of original reads
  nbKReads <- 0 #number of kept reads
  append=FALSE 
  allWin <- data.frame("Chr"=c(), "Start" = c(), "NbPositiveReads"= c(), "NbNegativeReads"= c(),"MaxCoverage" = c()) #data frame contains information of all windows for plotting
 
  for (chr in chromosomes){ #walk through each chromosome
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    mustKeepPos <- c()
    mustKeepNeg <- c()
    if (!missing(mustKeepRanges)){
      mkr <- mustKeepRanges[seqnames(mustKeepRanges)==chr]
      rg <- ranges(mkr)
      pos <- which(strand(mkr)=="+")
      mustKeepPos <- getStart(rg[pos,],win,step)  
      mustKeepNeg <- getStart(rg[-pos,],win,step)  
      remove(mkr)
      remove(rg)
      remove(pos)
    }
    alignmentInChr <- alignment[seqnames(alignment)==chr] #get the reads in the considering chromosome
    alignment <- alignment[seqnames(alignment)!=chr] #reduce the size of alignment, for memory purpose
    nbOReadsChr <- length(unique(names(alignmentInChr)))
    nbOReads <- nbOReads + nbOReadsChr
    
    results <- getKeptReadNames(alignmentInChr,covPos[[chromosomeIndex]],covNeg[[chromosomeIndex]],mustKeepPos,mustKeepNeg,getWin,readLength,len,win,step,pvalueThreshold,minCov,maxCov,logitThresholdP,logitThresholdM,errorRate)
    
    if (getWin){
      keptReadNames <- results$nameReads
      allWin <- rbind(allWin,results$Win %>% dplyr::mutate("Chr"=chr))
    }
    else{
      keptReadNames <- results
    }
    cat(paste0("Sequence ",chr,", length: ",len,", number of reads: ",nbOReadsChr,", number of kept reads: ",length(keptReadNames),"\n"),file=statfile,append=append)
    if (append==FALSE) append <- TRUE  
    if (length(keptReadNames)>0){##write the kept reads into output file
      #get the range of kept reads
      nbKReads <- nbKReads + length(keptReadNames)
      keptReads <- which(names(alignmentInChr) %in% keptReadNames)
      range <- rbamtools::bamRange(reader,c(chromosomeIndex-1,0,len))
      #write the kept reads into output file
      rbamtools::bamSave(writer,range[keptReads,],refid=chromosomeIndex-1)
      remove(range)
    }
    remove(keptReadNames)
  }
  rbamtools::bamClose(writer)
  remove(alignment)
  remove(covPos)
  remove(covNeg)
  rbamtools::bamClose(reader)
  
  cat("Summary:\n",file = statfile, append = append)
  cat(paste0("Number of original reads: ",nbOReads,", number of kept reads: ",nbKReads,", removal proportion: ",(nbOReads-nbKReads)/nbOReads),"\n",file = statfile,append=TRUE)
  endTime2 <- proc.time()
  cat(paste0("Total elapsed time ",(endTime2-startTime)[[3]]/60," minutes\n"),file = statfile,append=TRUE)
  if (getWin==TRUE){
    return(allWin)
  }
}
