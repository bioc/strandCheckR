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
  logitThresholdP <- binomial()$linkfun(threshold)
  logitThresholdM <- binomial()$linkfun(1-threshold)
  alignment <- GenomicAlignments::readGAlignments(bamfilein[1]) # read the first input alignment
  covPos <- alignment[strand(alignment)=="+"] %>% GenomicAlignments::coverage()
  covNeg <- alignment[strand(alignment)=="-"] %>% GenomicAlignments::coverage()
  
  lenSeq <- alignment@seqinfo@seqlengths #get the length of each chromosome as the number of bases
  allChromosomes <-alignment@seqinfo@seqnames #get the chromosome list
  if (missing(chromosomes)){
    chromosomes <- allChromosomes
  }
  emptyPosChromosomes <- sapply(covPos,function(ch){length(runValue(ch))==1 && runValue(ch)[1]==0})
  emptyNegChromosomes <- sapply(covNeg,function(ch){length(runValue(ch))==1 && runValue(ch)[1]==0})
  notEmptyChromosomes <- allChromosomes[!(emptyPosChromosomes * emptyNegChromosomes)]
  chromosomes <- intersect(chromosomes,notEmptyChromosomes)
  
  remove(alignment)
  if (length(bamfilein)>1){
    for (i in c(2:length(bamfilein))){
      alignment <- GenomicAlignments::readGAlignments(bamfilein[i])
      covPos <-covPos + alignment[strand(alignment)=="+"] %>% GenomicAlignments::coverage()
      covNeg <-covNeg + alignment[strand(alignment)=="-"] %>% GenomicAlignments::coverage()
      remove(alignment)
    }
  }
  keepWinPos <- list()
  keepWinNeg <- list()
  for (chr in chromosomes){
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    
    mustKeepPos <- c()
    mustKeepNeg <- c()
    if (!missing(mustKeepRanges)){
      mkr <- mustKeepRanges[seqnames(mustKeepRanges)==chr]
      rg <- ranges(mkr)
      pos <- which(strand(mkr)=="+")
      mustKeepPos <- getWin(rg[pos,],win,step)  
      mustKeepNeg <- getWin(rg[-pos,],win,step)  
      remove(mkr)
      remove(rg)
      remove(pos)
    }
    #compute the normalized value of each positive/negative window to be tested by pnorm
    windows <- computeWinInfo(runLength(covPos[[chromosomeIndex]]),runValue(covPos[[chromosomeIndex]]),runLength(covNeg[[chromosomeIndex]]),runValue(covNeg[[chromosomeIndex]]),readLength,len,win,step,minCov)
    
    plus <- dplyr::filter(windows,NbPositiveReads>=NbNegativeReads) %>% dplyr::mutate("propor" = NbPositiveReads/(NbPositiveReads+NbNegativeReads),"nbReads"=NbPositiveReads+NbNegativeReads)
    pvalueP <- pnorm(logitThresholdP,mean=binomial()$linkfun(plus$propor),sd=sqrt(1/(plus$nbReads)/plus$propor/(1-plus$propor)))
    
    minus <- dplyr::filter(windows,NbPositiveReads<NbNegativeReads) %>% dplyr::mutate("propor" = NbPositiveReads/(NbPositiveReads+NbNegativeReads),"nbReads"=NbPositiveReads+NbNegativeReads)
    pvalueM <- pnorm(logitThresholdM,mean=binomial()$linkfun(minus$propor),sd=sqrt(1/(minus$nbReads)/minus$propor/(1-minus$propor)),lower.tail = FALSE)
    
    if (maxCov>0){
      plus <- dplyr::mutate(plus,"pvalue"=pvalueP) %>% dplyr::filter(pvalueP<=pvalueThreshold | Start %in% mustKeepPos | MaxCoverage >= maxCov) %>% 
        dplyr::mutate("Start" = floor(Start/step)+1) %>%
        dplyr::select(c(Start,propor))
      minus <- dplyr::mutate(minus,"pvalue"=pvalueM) %>% dplyr::filter(pvalueM<=pvalueThreshold | Start %in% mustKeepNeg | MaxCoverage >= maxCov) %>% 
        dplyr::mutate("Start" = floor(Start/step)+1) %>%
        dplyr::select(c(Start,propor))
    }
    else{
      plus <- dplyr::mutate(plus,"pvalue"=pvalueP) %>% 
        dplyr::filter(pvalueP<=pvalueThreshold | Start %in% mustKeepPos) %>% 
        dplyr::mutate("Start" = floor(Start/step)+1) %>%
        dplyr::select(c(Start,propor))
      minus <- dplyr::mutate(minus,"pvalue"=pvalueM) %>% 
        dplyr::filter(pvalueM<=pvalueThreshold | Start %in% mustKeepNeg) %>% 
        dplyr::mutate("Start" = floor(Start/step)+1) %>%
        dplyr::select(c(Start,propor))
    }
    
    
    keepWinPos[[chromosomeIndex]] <- plus
    keepWinNeg[[chromosomeIndex]] <- minus
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
    alignment <- GenomicAlignments::readGAlignments(bamfilein[i],use.names = TRUE)
    reader <- rbamtools::bamReader(bamfilein[i],idx=TRUE)
    header <- rbamtools::getHeader(reader)
    writer <- rbamtools::bamWriter(header,bamfileout[i])
    for (chr in chromosomes){
      keptReads <- c()
      nbOReadsChr <- 0
      chromosomeIndex <- which(allChromosomes==chr)
      if (nrow(keepWinPos[[chromosomeIndex]])>0 || nrow(keepWinNeg[[chromosomeIndex]])>0){
        end <- lenSeq[chromosomeIndex]
        alignmentInChr <- alignment[seqnames(alignment)==chr]#get the reads in the considering chromosome
        alignment <- alignment[seqnames(alignment)!=chr]#reduce the size of alignment (for memory purpose)
        nbOReadsChr <- length(unique(names(alignmentInChr)))
        nbOReads[i] <- nbOReads[i] + nbOReadsChr
        
        strand <- as.vector(strand(alignmentInChr))
        indexPos <- which(strand=="+")
        indexNeg <- which(strand=="-")
        remove(strand)
        fragments <- getFragment(alignmentInChr)
        keptReads <- keepRead(fragments$Pos,fragments$Neg,keepWinPos[[chromosomeIndex]],keepWinNeg[[chromosomeIndex]],win,step,errorRate);   
        keptReads <- c(indexPos[unique(keptReads$Pos)],indexNeg[unique(keptReads$Neg)]) %>% sort() 
        keptReadNames <- names(alignmentInChr)[keptReads] %>% unique()
        keptReads <- which(names(alignmentInChr) %in% keptReadNames)
        remove(alignmentInChr)
        cat(paste0("Chromosome ",chr,", length: ",end,", number of original reads: ",nbOReadsChr,", number of kept reads: ",length(keptReadNames),"\n"),file=statfile,append=append)
        if (append==FALSE) append <- TRUE
      }
      if (length(keptReads)>0){
        nbKReads[i] <- nbKReads[i] + length(keptReadNames)
        range <- rbamtools::bamRange(reader,c(chromosomeIndex-1,0,end))
        rbamtools::bamSave(writer,range[keptReads,],refid=chromosomeIndex-1)
        remove(range)
        remove(keptReads)
      }
      if (i==length(bamfilein)){
        keepWinPos[[chromosomeIndex]] <- c(1)
        keepWinNeg[[chromosomeIndex]] <- c(1)
      }
    }
    remove(alignment)
    rbamtools::bamClose(reader)
    rbamtools::bamClose(writer)
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

