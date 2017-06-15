#' @title get the information of all windows from a paired end bam file
#'
#' @param bamfilein the input bam file\. Your bamfile should be sorted and have an index file located at the same path as well.
#' @param chromosomes the list of chromosomes to be read
#' @param yieldSize by default is 1e8, i.e. the bam file is read by block of chromosomes such that the total length of each block is at least 1e8
#' @param win the length of the sliding window, 1000 by default.
#' @param step the step length to sliding the window, 100 by default.
#' @param limit a read is considered to be included in a window if and only if at least limit percent of it is in the window. 0.75 by default
#' @param coverage if TRUE, then the strand information in each window corresponds to the sum of coverage coming from positive/negative reads; and not the number of positive/negative reads as default.
#'
#' @seealso filterDNA, filterDNAPairs, getWinFromBamFile, plotHist, plotWin
#' @export
#'
getWinFromPairedBamFile <- function(bamfilein,chromosomes,yieldSize=1e8,win=1000,step=100,limit=0.75,coverage=FALSE){
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
    bam <- scanBam(bamfilein,param=ScanBamParam(what=c("pos","cigar","strand","flag"),
                                                which=GRanges(seqnames = part,ranges = IRanges(start=1,end=lengthSeq[allChromosomes %in% part]))))
    nbOriginalReadsInChr <- sapply(1:length(bam),function(i){length(bam[[i]]$strand)})
    if (sum(nbOriginalReadsInChr)>0){
      nil <- which(nbOriginalReadsInChr!=0)
      nbOriginalReadsInChr <- nbOriginalReadsInChr[nil]
      part <- part[nil]
      nbOriginalReadsInPart <- c(0,cumsum(nbOriginalReadsInChr))
      lengthSeqInPart <- lengthSeqInPart[c(1,nil+1)]
      bam <- concatenateAlignments(bam[nil],nbOriginalReadsInPart,lengthSeqInPart,sum(nbOriginalReadsInChr),flag=TRUE)
      firstReadIndex <- ((floor(bam$flag/64) %% 2) == 1)
      secondReadIndex <- !firstReadIndex
      winFirstPositiveAlignments <- getWinOfAlignments(bam,"+",win,step,limit,firstReadIndex,coverage=coverage)
      winSecondPositiveAlignments <- getWinOfAlignments(bam,"+",win,step,limit,secondReadIndex,coverage=coverage)
      winFirstNegativeAlignments <- getWinOfAlignments(bam,"-",win,step,limit,firstReadIndex,coverage=coverage)
      winSecondNegativeAlignments <- getWinOfAlignments(bam,"-",win,step,limit,secondReadIndex,coverage=coverage)
      
      if (coverage){
        positiveFirstCoverage <- computeWinInfo(runLength(winFirstPositiveAlignments$Coverage),runValue(winFirstPositiveAlignments$Coverage),length(winFirstPositiveAlignments$Coverage),win,step)
        negativeFirstCoverage <- computeWinInfo(runLength(winFirstNegativeAlignments$Coverage),runValue(winFirstNegativeAlignments$Coverage),length(winFirstNegativeAlignments$Coverage),win,step)
        nbWinFirst <- max(max(positiveFirstCoverage$Start),max(negativeFirstCoverage$Start))
        nbFirstPositiveReads <- Rle(0,nbWin)
        nbFirstPositiveReads[positiveFirstCoverage$Start] <- positiveFirstCoverage$SumCoverage
        nbFirstNegativeReads <- Rle(0,nbWin)
        nbFirstNegativeReads[negativeFirstCoverage$Start] <- negativeFirstCoverage$SumCoverage
        
        positiveSecondCoverage <- computeWinInfo(runLength(winSecondPositiveAlignments$Coverage),runValue(winSecondPositiveAlignments$Coverage),length(winSecondPositiveAlignments$Coverage),win,step)
        negativeSecondCoverage <- computeWinInfo(runLength(winSecondNegativeAlignments$Coverage),runValue(winSecondNegativeAlignments$Coverage),length(winSecondNegativeAlignments$Coverage),win,step)
        nbWinSecond <- max(max(positiveSecondCoverage$Start),max(negativeSecondCoverage$Start))
        nbSecondPositiveReads <- Rle(0,nbWin)
        nbSecondPositiveReads[positiveSecondCoverage$Start] <- positiveSecondCoverage$SumCoverage
        nbSecondNegativeReads <- Rle(0,nbWin)
        nbSecondNegativeReads[negativeSecondCoverage$Start] <- negativeSecondCoverage$SumCoverage
      }
      else{
        nbFirstPositiveReads <- coverage(winFirstPositiveAlignments$Win)
        nbSecondPositiveReads <- coverage(winSecondPositiveAlignments$Win)
        nbFirstNegativeReads <- coverage(winFirstNegativeAlignments$Win)
        nbSecondNegativeReads <- coverage(winSecondNegativeAlignments$Win)
        lenFirstP <- length(nbFirstPositiveReads)
        lenSecondP <- length(nbSecondPositiveReads)
        lenFirstN <- length(nbFirstNegativeReads)
        lenSecondN <- length(nbSecondNegativeReads)
        if (lenFirstP>lenFirstN) {nbFirstNegativeReads <- c(nbFirstNegativeReads,rep(0,lenFirstP-lenFirstN))}
        else {nbFirstPositiveReads <- c(nbFirstPositiveReads,rep(0,lenFirstN-lenFirstP))}
        if (lenSecondP>lenSecondN) {nbSecondNegativeReads <- c(nbSecondNegativeReads,rep(0,lenSecondP-lenSecondN))}
        else {nbSecondPositiveReads <- c(nbSecondPositiveReads,rep(0,lenSecondN-lenSecondP))}
      }
      presentFirstWin <- which(as.vector((nbFirstPositiveReads>0) | (nbFirstNegativeReads>0))==TRUE)
      presentSecondWin <- which(as.vector((nbSecondPositiveReads>0) | (nbSecondNegativeReads>0))==TRUE)
      
      firstWin <- data.frame("Type"="First","Start" = presentFirstWin, "NbPositiveReads" = nbFirstPositiveReads[presentFirstWin], "NbNegativeReads" = nbFirstNegativeReads[presentFirstWin])
      secondWin <- data.frame("Type"="Second","Start" = presentSecondWin, "NbPositiveReads" = nbSecondPositiveReads[presentSecondWin], "NbNegativeReads" = nbSecondNegativeReads[presentSecondWin])
      ChromosomeFirst <- rep("",nrow(firstWin))
      ChromosomeSecond <- rep("",nrow(secondWin))
      for (i in seq_along(part)){
        mi <- ceiling((lengthSeqInPart[i]+1)/step)
        ma <- ceiling((lengthSeqInPart[i+1]-win+1)/step)
        
        j <- which(firstWin$Start >=mi & firstWin$Start <=ma)
        ChromosomeFirst[j] <- part[i]
        firstWin$Start[j] <- firstWin$Start[j] - mi +1
        
        j <- which(secondWin$Start >=mi & secondWin$Start <=ma)
        ChromosomeSecond[j] <- part[i]
        secondWin$Start[j] <- secondWin$Start[j] - mi +1
      }
      firstWin[["Chr"]] <- ChromosomeFirst
      secondWin[["Chr"]] <- ChromosomeSecond
      
      allWin <- rbind(allWin,firstWin)
      allWin <- rbind(allWin,secondWin)
    }
  }
  return(allWin)
}
