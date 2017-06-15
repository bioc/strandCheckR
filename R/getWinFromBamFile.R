#' @title get the information of all windows from a single end bam file
#'
#' @param bamfilein the input bam file\. Your bamfile should be sorted and have an index file located at the same path as well.
#' @param chromosomes the list of chromosomes to be read
#' @param yieldSize by default is 1e8, i.e. the bam file is read by block of chromosomes such that the total length of each block is at least 1e8
#' @param win the length of the sliding window, 1000 by default.
#' @param step the step length to sliding the window, 100 by default.
#' @param limit a read is considered to be included in a window if and only if at least limit percent of it is in the window. 0.75 by default
#' @param coverage if TRUE, then the strand information in each window corresponds to the sum of coverage coming from positive/negative reads; and not the number of positive/negative reads as default.
#' @seealso filterDNA, filterDNAPairs, getWinFromPairedBamFile, plotHist, plotWin
#' @export
#'
getWinFromBamFile <- function(bamfilein,chromosomes,yieldSize=1e8,win=1000,step=100,limit=0.75,coverage=FALSE){
  bf <- BamFile(bamfilein)
  seqinfo <- seqinfo(bf)
  allChromosomes <- seqnames(seqinfo)
  lengthSeq <- seqlengths(seqinfo)
  if (missing(chromosomes)){
    chromosomes <- allChromosomes
  }
  partition <- partitionChromosomes(chromosomes,lengthSeq[allChromosomes %in% chromosomes],yieldSize = yieldSize)
  allWin <- data.frame("Chr"=c(),"Start"=c(),"NbPositiveReads"=c(),"NbNegativeReads"=c())
  for (part in partition){
    lengthSeqInChr <- lengthSeq[allChromosomes %in% part]
    lengthSeqInPart <- c(0,cumsum(as.numeric(lengthSeqInChr)))
    lengthSeqInPart <- step*ceiling(lengthSeqInPart/step)
    
    bam <- scanBam(bamfilein,param=ScanBamParam(what=c("pos","cigar","strand"),
                                                which=GRanges(seqnames = part,ranges = IRanges(start=1,end=lengthSeqInChr))))
    
    nbOriginalReadsInChr <- sapply(1:length(bam),function(i){length(bam[[i]]$strand)})
    if (sum(nbOriginalReadsInChr)>0){
      nil <- which(nbOriginalReadsInChr!=0)
      nbOriginalReadsInChr <- nbOriginalReadsInChr[nil]
      part <- part[nil]
      nbOriginalReadsInPart <- c(0,cumsum(nbOriginalReadsInChr))
      lengthSeqInPart <- lengthSeqInPart[c(1,nil+1)]
      bam <- bamPart(bam[nil],nbOriginalReadsInPart,lengthSeqInPart,sum(nbOriginalReadsInChr))
      winPositiveAlignments <- getWinOfAlignments(bam,"+",win,step,limit,coverage=coverage)
      winNegativeAlignments <- getWinOfAlignments(bam,"-",win,step,limit,coverage=coverage)
      rm(bam)
      
      if (coverage){
        positiveCoverage <- computeWinInfo(runLength(winPositiveAlignments$Coverage),runValue(winPositiveAlignments$Coverage),length(winPositiveAlignments$Coverage),win,step)
        negativeCoverage <- computeWinInfo(runLength(winNegativeAlignments$Coverage),runValue(winNegativeAlignments$Coverage),length(winNegativeAlignments$Coverage),win,step)
        nbWin <- max(max(positiveCoverage$Start),max(negativeCoverage$Start))
        nbPositiveReads <- Rle(0,nbWin)
        nbPositiveReads[positiveCoverage$Start] <- positiveCoverage$SumCoverage
        nbNegativeReads <- Rle(0,nbWin)
        nbNegativeReads[negativeCoverage$Start] <- negativeCoverage$SumCoverage
      }
      else{
        nbPositiveReads <- coverage(winPositiveAlignments$Win)
        nbNegativeReads <- coverage(winNegativeAlignments$Win)
        lenP <- length(nbPositiveReads)
        lenN <- length(nbNegativeReads)
        if (lenP>lenN) {nbNegativeReads <- c(nbNegativeReads,rep(0,lenP-lenN))}
        else {nbPositiveReads <- c(nbPositiveReads,rep(0,lenN-lenP))}
      }
      
      presentWin <- which(as.vector((nbPositiveReads>0) | (nbNegativeReads>0))==TRUE)
      Win <- data.frame("Start" = presentWin, "NbPositiveReads" = nbPositiveReads[presentWin], "NbNegativeReads" = nbNegativeReads[presentWin])
      Chromosome <- rep("",nrow(Win))
      for (i in seq_along(part)){
        mi <- ceiling((lengthSeqInPart[i]+1)/step)
        ma <- ceiling((lengthSeqInPart[i+1]-win+1)/step)
        j <- which(Win$Start >=mi & Win$Start <=ma)
        Chromosome[j] <- part[i]
        Win$Start[j] <- Win$Start[j] - mi +1
      }
      Win[["Chr"]] <- Chromosome
      allWin <- rbind(allWin,Win)
    }
  }
  return(allWin %>% dplyr::mutate("Start" = (Start-1)*step+1))
}
