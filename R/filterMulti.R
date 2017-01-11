#' @title Filter Using Merged Files
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
#' 
#' @return the proportion of removed reads in every sample
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
  maxWin <- rep(0,length(idSeq))
  for (chr in chromosomes){
    chromosomeIndex <- which(idSeq==chr)
    len <- lenSeq[chromosomeIndex]
    #compute the normalized value of each positive/negative window to be tested by pnorm
    windows <- computeWinCov(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),readLength,len,win,step,minCov,maxCov,logitThreshold)
    
    maxWin[chromosomeIndex] <- 0
    if (nrow(windows$Plus)>0) maxWin[chromosomeIndex] <- windows$Plus$win[nrow(windows$Plus)]
    if (nrow(windows$Minus)>0) maxWin[chromosomeIndex] <- max(maxWin[chromosomeIndex],windows$Minus$win[nrow(windows$Minus)])
    
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
    alignment <- readGAlignments(bamfilein[i])
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
      
      position <- data.frame("group"=c(),"start"=c(),"end"=c(),"sample"=c()) #position of each fragment 
      alPos <- alignmentInChr[strand(alignmentInChr)=='+']
      alNeg <- alignmentInChr[strand(alignmentInChr)=='-']
      index <- getIndex(as.vector(strand(alignmentInChr)))
      remove(alignmentInChr)
      positionPos <- GenomicAlignments::extractAlignmentRangesOnReference(GenomicAlignments::cigar(alPos),pos=start(alPos)) %>% data.frame() %>% dplyr::select(-c(group_name,width))  
      positionPos <- positionPos[order(positionPos$start),]
      positionNeg <- GenomicAlignments::extractAlignmentRangesOnReference(GenomicAlignments::cigar(alNeg),pos=start(alNeg)) %>% data.frame() %>% dplyr::select(-c(group_name,width))  
      positionNeg <- positionNeg[order(positionNeg$start),]
      
      reads <- c()
      if (maxWin[chromosomeIndex]>0){
        reads <- keepRead(positionPos$start,positionPos$end,positionPos$group,positionNeg$start,positionNeg$end,positionNeg$group,keepWinPos[[chromosomeIndex]]$propor,keepWinNeg[[chromosomeIndex]]$propor,keepWinPos[[chromosomeIndex]]$win,keepWinNeg[[chromosomeIndex]]$win,maxWin[chromosomeIndex],win,step,errorRate);   
      }
      reads <- c(index$Pos[unique(reads$Pos)],index$Neg[unique(reads$Neg)]) %>% sort() 
      if (chromosomeIndex>1){
        append <- TRUE
      }
      cat(paste0("Chromosome ",chr,", length: ",end,", number of original reads: ",nbOReadsChr,", number of kept reads: ",length(reads),"\n"),file=statfile,append=append)
      nbKReads[i] <- nbKReads[i] + length(reads)
      if (length(reads)>0){
        range <- bamRange(reader,c(chromosomeIndex-1,0,end))
        bamSave(writer,range[reads,],refid=chromosomeIndex-1)
        remove(range)
      }
      remove(reads)
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

