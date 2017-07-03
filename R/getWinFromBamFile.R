#' @title get the number of positive/negative reads of all windows from a single end bam file
#'
#' @param bamfilein the input single-end bam file. Your bamfile should be sorted and have an index file located at the same path as well.
#' @param chromosomes the list of chromosomes to be read
#' @param yieldSize by default is 1e8, i.e. the bam file is read by block of chromosomes such that the total length of each block is at least 1e8
#' @param win the length of the sliding window, 1000 by default.
#' @param step the step length to sliding the window, 100 by default.
#' @param limit a read is considered to be included in a window if and only if at least \code{limit} percent of it is in the window. 0.75 by default
#' @param coverage if TRUE, then the strand information in each window corresponds to the sum of coverage coming from positive/negative reads; and not the number of positive/negative reads as default.
#' @seealso filterDNA, filterDNAPairs, getWinFromPairedBamFile, plotHist, plotWin
#' @export
#' @importFrom IRanges Views
#' @examples 
#' bamfilein <- system.file("data","s1.chr1.bam",package = "rnaCleanR")
#' win <- getWin(bamfilein)
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
  if (coverage==TRUE){
    allWin <- data.frame("Chr"=c(), "Start" = c(), "NbPositiveReads"= c(), "NbNegativeReads"= c(),"MaxCoverage" = c()) 
  } else{
    allWin <- data.frame("Chr"=c(), "Start" = c(), "NbPositiveReads"= c(), "NbNegativeReads"= c())  
  }
  statInfo <- data.frame("Sequence"="chr","Length"=rep(0,length(chromosomes)),
                         "NbOriginalReads" = rep(0,length(chromosomes)), 
                         "NbKeptReads" = rep(0,length(chromosomes)),
                         "FirstBaseInPartition" = rep(NA,length(chromosomes)),
                         "LastBaseInPartition" = rep(NA,length(chromosomes)),
                         "FirstReadInPartition" = rep(NA,length(chromosomes)),
                         "LastReadInPartition" = rep(NA,length(chromosomes)),
                         stringsAsFactors = FALSE)
  
  for (part in partition){
    idPart <- which(chromosomes %in% part)
    statInfo$Sequence[idPart] <- part
    statInfo$Length[idPart] <- lengthSeq[allChromosomes %in% part]
    
    bam <- scanBam(bamfilein,param=ScanBamParam(what=c("pos","cigar","strand"),
                                                which=GRanges(seqnames = part,ranges = IRanges(start=1,end=statInfo$Length[idPart]))))
    statInfo$NbOriginalReads[idPart] <- sapply(seq_along(bam),function(i){length(bam[[i]]$strand)})
    if (sum(statInfo$NbOriginalReads[idPart])>0){
      statInfo[idPart,] <- statInfoInPartition(statInfo[idPart,],step)
      bam <- concatenateAlignments(bam,statInfo[idPart,])
      winPositiveAlignments <- getWinOfAlignments(bam,"+",win,step,limit,coverage=coverage)
      winNegativeAlignments <- getWinOfAlignments(bam,"-",win,step,limit,coverage=coverage)
      rm(bam)
      if (coverage){
        lenPC <- length(winPositiveAlignments$Coverage)
        lenNC <- length(winNegativeAlignments$Coverage)
        if (lenPC>lenNC) {winNegativeAlignments$Coverage <- c(winNegativeAlignments$Coverage,rep(0,lenPC-lenNC))} else {winPositiveAlignments$Coverage <- c(winPositiveAlignments$Coverage,rep(0,lenNC-lenPC))}
        nbWin <- ceiling((length(winPositiveAlignments$Coverage)-win)/step)+1
        nbPositiveReads <- Views(winPositiveAlignments$Coverage,
                                 start = seq(1,(nbWin-1)*step+1,step),
                                 end=seq(win,(nbWin-1)*step+win,step)) %>%
          sum() %>% Rle()
        nbNegativeReads <- Views(winNegativeAlignments$Coverage,
                                 start = seq(1,(nbWin-1)*step+1,step),
                                 end=seq(win,(nbWin-1)*step+win,step)) %>%
          sum() %>% Rle()
        maxCoverage <- Views(winPositiveAlignments$Coverage+winNegativeAlignments$Coverage,
                             start = seq(1,(nbWin-1)*step+1,step),
                             end=seq(win,(nbWin-1)*step+win,step)) %>% 
          max() %>% Rle()
      }
      else{
        nbPositiveReads <- coverage(winPositiveAlignments$Win)
        nbNegativeReads <- coverage(winNegativeAlignments$Win)
        lenP <- length(nbPositiveReads)
        lenN <- length(nbNegativeReads)
        if (lenP>lenN) {nbNegativeReads <- c(nbNegativeReads,rep(0,lenP-lenN))} else {nbPositiveReads <- c(nbPositiveReads,rep(0,lenN-lenP))}
      }
      presentWin <- which(as.vector((nbPositiveReads>0) | (nbNegativeReads>0))==TRUE)
      if (coverage){
        Win <- data.frame("Start" = presentWin, "NbPositiveReads" = nbPositiveReads[presentWin], "NbNegativeReads" = nbNegativeReads[presentWin],"MaxCoverage" = maxCoverage[presentWin])  
      }
      else{
        Win <- data.frame("Start" = presentWin, "NbPositiveReads" = nbPositiveReads[presentWin], "NbNegativeReads" = nbNegativeReads[presentWin])  
      }
      Chromosome <- rep("",nrow(Win))
      for (i in seq_along(part)){
        id <- which(chromosomes == part[i])
        mi <- ceiling(statInfo$FirstBaseInPartition[id]/step)
        ma <- ceiling((statInfo$LastBaseInPartition[id]-win+1)/step)
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
