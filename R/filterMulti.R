#' @title Filter Several Bam Files
#' 
#' @description Filter several bamfiles based on the strand specific of the merged files. 
#' The strand of each sliding window is caculated based on coverage.
#' 
#' @details None
#' 
#' @param bamfilein the input bam files to be filterd
#' @param bamfileout the output bam files
#' @param chromosomes the list of chromosomes to be filtered
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' @param threshold the threshold upper which we keep the reads
#' @param pvalueThreshold the threshold for the p-value
#' @param minR if a window has the max coverage least than minR, then it will be rejected
#' @param maxR if a window has the max coverage greater than maxR, then it will be kept
#'
#' @details filterMulti reads a set of bam files containing strand specific RNA reads, and filter the potential double strand contamination DNA from all these files. 
#' This method also uses sliding windows approach as method filterOne, but it uses the strand information of the reads coming from all input bam files.
#' 
#' @examples  
#' bamfilein <- system.file("data",c("s1.chr1.bam","s2.chr1.bam"),package = "rnaCleanR")
#' filterMulti(bamfilein,bamfileout=c("s1.filter.bam","s2.filter.bam"),statfile = "out.stat",histPlot = TRUE,histfile = "hist.pdf", readLength = 100, threshold = 0.7)
#' @export
#' 
filterMulti <- function(bamfilein,bamfileout,statfile,chromosomes=NULL,readLength,win,step,pvalueThreshold=0.05,minCov=0,maxCov=0,threshold=0.7,errorRate=0.01){
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
    message("Summary is written into file out.stat")
    statfile <- "out.stat"
  }
  logitThreshold <- binomial()$linkfun(threshold)
  alignment <- GenomicAlignments::readGAlignments(bamfilein[1]) # read the first input alignment
  covPos <- alignment[strand(alignment)=="+"] %>% GenomicAlignments::coverage()
  covNeg <- alignment[strand(alignment)=="-"] %>% GenomicAlignments::coverage()
  idSeq <- levels(seqnames(alignment)) #get the chromosome list
  remove(alignment)
  if (is.null(chromosomes)) chromosomes <- idSeq
  lenSeq<-sapply(1:length(covPos),function(i) length(covPos[[i]])) #get the length of each chromosome as the number of bases
  if (length(bamfilein)>1){
    for (i in c(2:length(bamfilein))){
      alignment <- readGAlignments(bamfilein[i])
      covPos <-covPos +alignment[strand(alignment)=="+"] %>% GenomicAlignments::coverage()
      covNeg <-covNeg+alignment[strand(alignment)=="-"] %>% GenomicAlignments::coverage()
      remove(alignment)
    }
  }
  keepWinPos <- list()
  keepWinNeg <- list()
  for (chr in chromosomes){
    chromosomeIndex <- which(idSeq==chr)
    len <- lenSeq[chromosomeIndex]
    #compute the normalized value of each positive/negative window to be tested by pnorm
    windows <- computeWin(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),readLength,len,win,step,minCov,maxCov,logitThreshold)
    
    windows$Plus <- dplyr::mutate(windows$Plus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% dplyr::filter( pvalue <= pvalueThreshold) %>% dplyr::select(c(win,propor))
    windows$Minus <- dplyr::mutate(windows$Minus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% dplyr::filter(pvalue <= pvalueThreshold) %>% dplyr::select(c(win,propor))
    keepWinPos[[chromosomeIndex]] <- windows$Plus
    keepWinNeg[[chromosomeIndex]] <- windows$Minus
    remove(windows)
  }
  remove(covPos)
  remove(covNeg)
  
  #write the kept reads to output files
  nbOReads <- rep(0,length(bamfilein))
  nbKReads <- rep(0,length(bamfilein))
  append <- FALSE
  for (i in seq_along(bamfilein)){
    cat(paste0("File ",i," : ",bamfilein[i],"\n"),file = statfile,append = append)
    alignment <- GenomicAlignments::readGAlignments(bamfilein[i])
    reader <- bamReader(bamfilein[i],idx=TRUE)
    header <- getHeader(reader)
    writer <- bamWriter(header,bamfileout[i])
    if (i>1){
      append <- TRUE
    }
    for (chr in chromosomes){
      chromosomeIndex <- which(idSeq==chr)
      end <- lenSeq[chromosomeIndex]
      alignmentInChr <- alignment[seqnames(alignment)==chr]#get the reads in the considering chromosome
      alignment <- alignment[seqnames(alignment)!=chr]#reduce the size of alignment (for memory purpose)
      nbOReadsChr <- length(alignmentInChr)
      nbOReads[i] <- nbOReads[i] + nbOReadsChr
      
      index <- getIndex(as.vector(strand(alignmentInChr)))
      fragments <- getFragment(alignmentInChr)
      remove(alignmentInChr)
      
      keptReads <- c()
      if (nrow(keepWinPos[[chromosomeIndex]])>0 || nrow(keepWinNeg[[chromosomeIndex]])>0){
        keptReads <- keepRead(fragments$Pos,fragments$Neg,keepWinPos[[chromosomeIndex]],keepWinNeg[[chromosomeIndex]],win,step,errorRate);   
      }
      keptReads <- c(index$Pos[unique(keptReads$Pos)],index$Neg[unique(keptReads$Neg)]) %>% sort() 
      remove(keptFragments)
      if (chromosomeIndex>1){
        append <- TRUE
      }
      cat(paste0("Chromosome ",chr,", length: ",end,", number of original reads: ",nbOReadsChr,", number of kept reads: ",length(keptReads),"\n"),file=statfile,append=append)
      nbKReads[i] <- nbKReads[i] + length(keptReads)
      if (length(keptReads)>0){
        range <- bamRange(reader,c(chromosomeIndex-1,0,end))
        bamSave(writer,range[keptReads,],refid=chromosomeIndex-1)
        remove(range)
      }
      remove(keptReads)
      if (i==length(bamfilein)){
        keepWinPos[[chromosomeIndex]] <- c(1)
        keepWinNeg[[chromosomeIndex]] <- c(1)
      }
    }
    remove(alignment)
    bamClose(reader)
    bamClose(writer)
  }
  remove(keepWinPos)
  remove(keepWinNeg)
  endTime2 <- proc.time()
  cat("Summary:\n",file = statfile, append = append)
  for (i in seq_along(bamfilein)){
    cat(paste0("File ",i," ",bamfilein[i]," - Number of original reads: ",nbOReads[i],". Number of kept reads: ",nbKReads[i],". Removal proportion: ",(nbOReads[i]-nbKReads[i])/nbOReads[i]),"\n",file = statfile,append = append)
  }
  cat(paste0("Total elapsed time ",(endTime2-startTime)[[3]]/60," minutes"),file = statfile,append = append)
}

