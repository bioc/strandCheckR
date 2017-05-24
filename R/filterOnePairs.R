#' @title Filter One Paried End Bam File 
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
#' filterOne(bamfilein,bamfileout="out.bam",statfile = "out.stat", readLength = 100, threshold = 0.7)
#' @export
#' 
filterOnePairs <- function(bamfilein,bamfileout,statfile,chromosomes,mustKeepRanges,getWin = FALSE,readLength,win,step,pvalueThreshold=0.05,minCov=0,maxCov=0,threshold=0.7,errorRate=0.01,pair="free"){
  startTime <- proc.time()
  if (!(pair %in% c("free","intersect","union"))){
    stop("pair should be either free, or intersect, or union")
  }
  
  # read the input alignments and compute positive/negative coverge
  alignment <- GenomicAlignments::readGAlignments(bamfilein,use.names = TRUE,param = ScanBamParam(what="flag"))
  firstReadIndex <- ((floor(alignment@elementMetadata$flag/64) %% 2) == 1) 
  alignmentFirst <- alignment[firstReadIndex]
  alignmentSecond <- alignment[!firstReadIndex]  

  
  covPosFirst <- alignmentFirst[strand(alignmentFirst)=="+"] %>% GenomicAlignments::coverage() 
  covNegFirst <- alignmentFirst[strand(alignmentFirst)=="-"] %>% GenomicAlignments::coverage()
  
  covPosSecond <- alignmentSecond[strand(alignmentSecond)=="+"] %>% GenomicAlignments::coverage() 
  covNegSecond <- alignmentSecond[strand(alignmentSecond)=="-"] %>% GenomicAlignments::coverage() 
  
  rm(firstReadIndex)
  rm(alignmentFirst)
  rm(alignmentSecond)
  
  allChromosomes <- alignment@seqinfo@seqnames #get the name of each chromosome 
  lenSeq <- alignment@seqinfo@seqlengths #get the length of each chromosome 
  
  
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
  
  if (missing(chromosomes)){
    chromosomes <- allChromosomes
  }  
  emptyPosChromosomes1 <- sapply(covPosFirst,function(ch){length(runValue(ch))==1 && runValue(ch)[1]==0})
  emptyPosChromosomes2 <- sapply(covPosSecond,function(ch){length(runValue(ch))==1 && runValue(ch)[1]==0})
  emptyNegChromosomes1 <- sapply(covNegFirst,function(ch){length(runValue(ch))==1 && runValue(ch)[1]==0})
  emptyNegChromosomes2 <- sapply(covNegSecond,function(ch){length(runValue(ch))==1 && runValue(ch)[1]==0})
  notEmptyChromosomes <- allChromosomes[!(emptyPosChromosomes1 * emptyPosChromosomes2 * emptyNegChromosomes1 * emptyNegChromosomes2)]
  
  chromosomes <- intersect(chromosomes,notEmptyChromosomes)
  
  reader <- rbamtools::bamReader(bamfilein,idx=TRUE) #open a reader of the input bamfile to extract read afterward
  header <- rbamtools::getHeader(reader) #get the header of the input bam file
  writer <- rbamtools::bamWriter(header,bamfileout) #prepare to write the output bamfile with the same header
  logitThresholdP <- binomial()$linkfun(threshold)
  logitThresholdM <- binomial()$linkfun(1-threshold)
  nbOFirstReads <- 0 #number of original reads
  nbKFirstReads <- 0 #number of kept reads
  nbOSecondReads <- 0 #number of original reads
  nbKSecondReads <- 0 #number of kept reads
  append=FALSE 
  allWinFirst <- data.frame("Chr"=c(),"Proportion"=c(),"NbReads"=c(),"MaxCoverage"=c()) 
  allWinSecond <- data.frame("Chr"=c(),"Proportion"=c(),"NbReads"=c(),"MaxCoverage"=c())
  
  for (chr in chromosomes){ #walk through each chromosome
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    
    mustKeepPosWin <- c()
    mustKeepNegWin <- c()
    if (!missing(mustKeepRanges)){
      mkr <- mustKeepRanges[seqnames(mustKeepRanges)==chr] 
      pos <- which(strand(mkr)=="+")
      rg <- ranges(mkr)
      mustKeepPosWin <- getWin(rg[pos,],win,step)  
      mustKeepNegWin <- getWin(rg[-pos,],win,step) 
      remove(mkr)
      remove(rg)
      remove(pos)
    }
    
    alignmentInChr <- alignment[seqnames(alignment)==chr]
    alignment <- alignment[seqnames(alignment)!=chr]
    firstReadIndex <- ((floor(alignmentInChr@elementMetadata$flag/64) %% 2) == 1) 
    alignmentInChrFirst <- alignmentInChr[firstReadIndex]
    alignmentInChrSecond <- alignmentInChr[!firstReadIndex]  
    
    nbOFirstReadsChr <- length(unique(names(alignmentInChrFirst)))
    nbOFirstReads <- nbOFirstReads + nbOFirstReadsChr
    
    nbOSecondReadsChr <- length(unique(names(alignmentInChrSecond)))
    nbOSecondReads <- nbOSecondReads + nbOSecondReadsChr
    
    first <- getKeptReadNames(alignmentInChrFirst,covPosFirst[[chromosomeIndex]],covNegFirst[[chromosomeIndex]],mustKeepPosWin,mustKeepNegWin,getWin,readLength,len,win,step,pvalueThreshold,minCov,maxCov,logitThresholdP,logitThresholdM,errorRate)
    
    second <- getKeptReadNames(alignmentInChrSecond,covPosSecond[[chromosomeIndex]],covNegSecond[[chromosomeIndex]],mustKeepPosWin,mustKeepNegWin,getWin,readLength,len,win,step,pvalueThreshold,minCov,maxCov,logitThresholdP,logitThresholdM,errorRate)
  
    if (getWin){
      keptFirstReadNames <- first$nameReads
      keptSecondReadNames <- second$nameReads
      allWinFirst <- rbind(allWinFirst,first$Win %>% dplyr::mutate("Type"="First","Chr"=chr))
      allWinSecond <- rbind(allWinSecond,second$Win %>% dplyr::mutate("Type"="Second","Chr"=chr))
    }
    else{
      keptFirstReadNames <- first
      keptSecondReadNames <- second
    }
    rm(first)
    rm(second)
    if (pair=="free"){
      nbKFirstReadsChr <- length(keptFirstReadNames)
      nbKSecondReadsChr <- length(keptSecondReadNames)
      cat(paste0("Sequence ",chr,", length: ",len,", number of first reads: ",nbOFirstReadsChr,", number of second reads: ",nbOSecondReadsChr,", number of kept first reads: ",nbKFirstReadsChr,", number of kept second reads: ",nbKSecondReadsChr,"\n"),file=statfile,append=append)    
      keptFirstAlignments <- which(firstReadIndex==TRUE)[which(names(alignmentInChrFirst) %in% keptFirstReadNames)]
      keptSecondAlignments <- which(firstReadIndex==FALSE)[which(names(alignmentInChrSecond) %in% keptSecondReadNames)]
      nbKFirstReads <- nbKFirstReads + nbKFirstReadsChr
      nbKSecondReads <- nbKSecondReads + nbKSecondReadsChr
      keptReads <-  c(keptFirstAlignments,keptSecondAlignments) %>% sort()
    }
    else{
      if (pair=="intersect"){
        keptReadNames <- intersect(keptFirstReadNames,keptSecondReadNames)
      }
      if (pair=="union"){
        commonNames <- intersect(names(alignmentInChrFirst),names(alignmentInChrSecond))
        keptReadNames <- intersect(commonNames,union(keptFirstReadNames,keptSecondReadNames))
      }
      nbKFirstReads <- nbKFirstReads + length(keptReadNames)
      nbKSecondReads <- nbKSecondReads + length(keptReadNames)
      cat(paste0("Sequence ",chr,", length: ",len,", number of first reads: ",nbOFirstReadsChr,", number of second reads: ",nbOSecondReadsChr,", number of kept paired reads: ",length(keptReadNames),"\n"),file=statfile,append=append)    
      keptReads <- which(names(alignmentInChr) %in% keptReadNames)
    } 
    
    if (append==FALSE) append <- TRUE  
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
  rbamtools::bamClose(reader)
  
  cat("Summary:\n",file = statfile, append = append)
  cat(paste0("Number of original first reads: ",nbOFirstReads,", number of original second reads: ",nbOSecondReads,", number of kept first reads: ",nbKFirstReads,", number of kept second reads: ",nbKSecondReads,", removal proportion of first reads: ",(nbOFirstReads-nbKFirstReads)/nbOFirstReads),", removal proportion of second reads: ",(nbOSecondReads-nbKSecondReads)/nbOSecondReads,"\n",file = statfile,append=TRUE)
  endTime2 <- proc.time()
  cat(paste0("Total elapsed time ",(endTime2-startTime)[[3]]/60," minutes\n"),file = statfile,append=TRUE)
  if (getWin){
    return(rbind(allWinFirst,allWinSecond))
  }
}
