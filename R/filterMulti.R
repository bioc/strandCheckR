#' @title Filter Several Bam Files
#' 
#' @description Filter several bam files based on the strand specific of the merged files. 
#' The strand of each sliding window is calculated based on coverage.
#' 
#' 
#' @param bamfilein the input bam files to be filterd
#' @param bamfileout the output bam files
#' @param chromosomes the list of chromosomes to be filtered
#' @param mustKeepRanges a GRanges object defines the ranges such that every read maps to those ranges must be always kept regardless the strand proportion of the windows containing them.
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' @param threshold the threshold upper which we keep the reads
#' @param pvalueThreshold the threshold for the p-value
#' @param minR if a window has the max coverage least than minR, then it will be rejected
#' @param maxR if a window has the max coverage greater than maxR, then it will be kept
#'
#' @details filterMulti reads a set of bam files containing strand specific RNA reads, and filter the potential double strand contamination DNA from all these files. 
#' Similar to the filterOne method, it uses a window sliding across the whole genome, but on the merged alignments from all input bam file.
#' See filterOne for further information.
#' 
#' @seealso filterOne, getPlot
#' 
#' @examples  
#' bamfileins <- system.file("data",c("s1.chr1.bam","s2.chr1.bam"),package = "rnaCleanR")
#' filterMulti(bamfileins,bamfileout=c("s1.filter.bam","s2.filter.bam"),statfile = "out.stat", readLength = 100, threshold = 0.7)
#' @export
#' 
filterMulti <- function(bamfilein,bamfileout,statfile,chromosomes,mustKeepRanges,readLength,win,step,pvalueThreshold=0.05,minCov=0,maxCov=0,threshold=0.7,errorRate=0.01){
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
  lenSeq <- alignment@seqinfo@seqlengths #get the length of each chromosome as the number of bases
  allChromosomes <-alignment@seqinfo@seqnames #get the chromosome list
  
  covPos <- alignment[strand(alignment)=="+"] %>% GenomicAlignments::coverage()
  covNeg <- alignment[strand(alignment)=="-"] %>% GenomicAlignments::coverage()
  
  remove(alignment)
  if (missing(chromosomes)){
    chromosomes <- allChromosomes
  }  
  
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
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    #compute the normalized value of each positive/negative window to be tested by pnorm
    windows <- computeWin(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),readLength,len,win,step,minCov,maxCov,logitThreshold)
    if (!missing(mustKeepRanges)){
      mkr <- mustKeepRanges[seqnames(mustKeepRanges)==chr]
      rg <- ranges(mkr)
      pos <- which(strand(mkr)=="+")
      mustKeepPosWin <- getWin(rg[pos,],win,step)  
      mustKeepNegWin <- getWin(rg[-pos,],win,step)  
      windows$Plus <- dplyr::mutate(windows$Plus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% dplyr::filter( pvalue <= pvalueThreshold | win %in% mustKeepPosWin) %>% dplyr::select(c(win,propor))
      windows$Minus <- dplyr::mutate(windows$Minus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% dplyr::filter(pvalue <= pvalueThreshold | win %in% mustKeepNegWin) %>% dplyr::select(c(win,propor))
      remove(mkr)
      remove(rg)
      remove(pos)
    }
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
      keptReads <- c()
      nbOReadsChr <- 0
      chromosomeIndex <- which(allChromosomes==chr)
      if (nrow(keepWinPos[[chromosomeIndex]])>0 || nrow(keepWinNeg[[chromosomeIndex]])>0){
        chromosomeIndex <- which(allChromosomes==chr)
        end <- lenSeq[chromosomeIndex]
        alignmentInChr <- alignment[seqnames(alignment)==chr]#get the reads in the considering chromosome
        alignment <- alignment[seqnames(alignment)!=chr]#reduce the size of alignment (for memory purpose)
        nbOReadsChr <- length(alignmentInChr)
        nbOReads[i] <- nbOReads[i] + nbOReadsChr
        
        strand <- as.vector(strand(alignmentInChr))
        indexPos <- which(strand=="+")
        indexNeg <- which(strand=="-")
        remove(strand)
        fragments <- getFragment(alignmentInChr)
        remove(alignmentInChr)
        keptReads <- keepRead(fragments$Pos,fragments$Neg,keepWinPos[[chromosomeIndex]],keepWinNeg[[chromosomeIndex]],win,step,errorRate);   
        keptReads <- c(indexPos[unique(keptReads$Pos)],indexNeg[unique(keptReads$Neg)]) %>% sort() 
        
        cat(paste0("Chromosome ",chr,", length: ",end,", number of original reads: ",nbOReadsChr,", number of kept reads: ",length(keptReads),"\n"),file=statfile,append=append)
        if (append==FALSE) append <- TRUE
      }
      if (length(keptReads)>0){
        nbKReads[i] <- nbKReads[i] + length(keptReads)
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

