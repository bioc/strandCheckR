#' @title get the number of positive/negatives of all windows from a paired end bam file
#'
#' @param bamfilein the input paired-end bam file. Your bamfile should be sorted and have an index file located at the same path as well.
#' @param chromosomes the list of chromosomes to be read
#' @param yieldSize by default is 1e8, i.e. the bam file is read by block of chromosomes such that the total length of each block is at least 1e8
#' @param win the length of the sliding window, 1000 by default.
#' @param step the step length to sliding the window, 100 by default.
#' @param limit a read is considered to be included in a window if and only if at least \code{limit} percent of it is in the window. 0.75 by default
#' @param coverage if TRUE, then the strand information in each window corresponds to the sum of coverage coming from positive/negative reads; and not the number of positive/negative reads as default.
#'
#' @seealso filterDNA, filterDNAPairs, getWinFromBamFile, plotHist, plotWin
#' @export
#' @examples
#' bamfilein <- system.file("data","120.10.bam",package = "rnaCleanR")
#' win <- getWinFromPairedBamFile(bamfilein)
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
  if (coverage==TRUE){
    allWin <- data.frame("Type"=c(),"Chr"=c(), "Start" = c(), "NbPositiveReads"= c(), "NbNegativeReads"= c(),"MaxCoverage" = c()) 
  }
  else{
    allWin <- data.frame("Type"=c(),"Chr"=c(), "Start" = c(), "NbPositiveReads"= c(), "NbNegativeReads"= c())  
  }
  statInfo <- data.frame("Sequence"="chr","Length"=rep(0,length(chromosomes)),
                         "NbOriginalReads" = rep(0,length(chromosomes)), 
                         "NbOriginalFirstReads" = rep(0,length(chromosomes)), 
                         "NbOriginalSecondReads" = rep(0,length(chromosomes)), 
                         "NbKeptFirstReads" = rep(0,length(chromosomes)),
                         "NbKeptSecondReads" = rep(0,length(chromosomes)),
                         "FirstBaseInPartition" = rep(NA,length(chromosomes)),
                         "LastBaseInPartition" = rep(NA,length(chromosomes)),
                         "FirstReadInPartition" = rep(NA,length(chromosomes)),
                         "LastReadInPartition" = rep(NA,length(chromosomes)),
                         stringsAsFactors = FALSE)
  
  for (part in partition){
    idPart <- which(chromosomes %in% part)
    statInfo$Sequence[idPart] <- part
    statInfo$Length[idPart] <- lengthSeq[allChromosomes %in% part]
    
    bam <- scanBam(bamfilein,param=ScanBamParam(what=c("pos","cigar","strand","flag"),
                                                which=GRanges(seqnames = part,ranges = IRanges(start=1,end=statInfo$Length[idPart]))))
    statInfo$NbOriginalReads[idPart] <- sapply(seq_along(bam),function(i){length(bam[[i]]$strand)})
    if (sum(statInfo$NbOriginalReads[idPart])>0){
      statInfo[idPart,] <- statInfoInPartition(statInfo[idPart,],step)
      bam <- concatenateAlignments(bam,statInfo[idPart,],flag=TRUE)
      
      firstReadIndex <- ((floor(bam$flag/64) %% 2) == 1)
      secondReadIndex <- !firstReadIndex
      winFirstPositiveAlignments <- getWinOfAlignments(bam,"+",win,step,limit,firstReadIndex,coverage=coverage)
      winSecondPositiveAlignments <- getWinOfAlignments(bam,"+",win,step,limit,secondReadIndex,coverage=coverage)
      winFirstNegativeAlignments <- getWinOfAlignments(bam,"-",win,step,limit,firstReadIndex,coverage=coverage)
      winSecondNegativeAlignments <- getWinOfAlignments(bam,"-",win,step,limit,secondReadIndex,coverage=coverage)
      
      if (coverage){
        lenPC <- length(winFirstPositiveAlignments$Coverage)
        lenNC <- length(winFirstNegativeAlignments$Coverage)
        if (lenPC>lenNC) {winFirstNegativeAlignments$Coverage <- c(winFirstNegativeAlignments$Coverage,rep(0,lenPC-lenNC))} else {winFirstPositiveAlignments$Coverage <- c(winFirstPositiveAlignments$Coverage,rep(0,lenNC-lenPC))}
        nbWin <- ceiling((length(winFirstPositiveAlignments$Coverage)-win)/step)+1
        nbFirstPositiveReads <- Views(winFirstPositiveAlignments$Coverage,
                                      start = seq(1,(nbWin-1)*step+1,step),
                                      end=seq(win,(nbWin-1)*step+win,step)) %>%
          sum() %>% Rle()
        nbFirstNegativeReads <- Views(winFirstNegativeAlignments$Coverage,
                                      start = seq(1,(nbWin-1)*step+1,step),
                                      end=seq(win,(nbWin-1)*step+win,step)) %>%
          sum() %>% Rle()
        maxFirstCoverage <- Views(winFirstPositiveAlignments$Coverage+winFirstNegativeAlignments$Coverage,
                                  start = seq(1,(nbWin-1)*step+1,step),
                                  end=seq(win,(nbWin-1)*step+win,step)) %>%
          max() %>% Rle()
        
        lenPC <- length(winSecondPositiveAlignments$Coverage)
        lenNC <- length(winSecondNegativeAlignments$Coverage)
        nbWin <- ceiling((length(winSecondPositiveAlignments$Coverage)-win)/step)+1
        if (lenPC>lenNC) {winSecondNegativeAlignments$Coverage <- c(winSecondNegativeAlignments$Coverage,rep(0,lenPC-lenNC))} else {winSecondPositiveAlignments$Coverage <- c(winSecondPositiveAlignments$Coverage,rep(0,lenNC-lenPC))}
        nbSecondPositiveReads <- Views(winSecondPositiveAlignments$Coverage,
                                       start = seq(1,(nbWin-1)*step+1,step),
                                       end=seq(win,(nbWin-1)*step+win,step)) %>%
          sum() %>% Rle()
        nbSecondNegativeReads <- Views(winSecondNegativeAlignments$Coverage,
                                       start = seq(1,(nbWin-1)*step+1,step),
                                       end=seq(win,(nbWin-1)*step+win,step)) %>%
          sum() %>% Rle()
        maxSecondCoverage <- Views(winSecondPositiveAlignments$Coverage+winSecondNegativeAlignments$Coverage,
                                   start = seq(1,(nbWin-1)*step+1,step),
                                   end=seq(win,(nbWin-1)*step+win,step)) %>%
          max() %>% Rle()
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
        if (lenFirstP>lenFirstN) {
          nbFirstNegativeReads <- c(nbFirstNegativeReads,rep(0,lenFirstP-lenFirstN))
        } else {
          nbFirstPositiveReads <- c(nbFirstPositiveReads,rep(0,lenFirstN-lenFirstP))
        }
        if (lenSecondP>lenSecondN) {
          nbSecondNegativeReads <- c(nbSecondNegativeReads,rep(0,lenSecondP-lenSecondN))
        } else {
            nbSecondPositiveReads <- c(nbSecondPositiveReads,rep(0,lenSecondN-lenSecondP))
        }
      }
      presentFirstWin <- which(as.vector((nbFirstPositiveReads>0) | (nbFirstNegativeReads>0))==TRUE)
      presentSecondWin <- which(as.vector((nbSecondPositiveReads>0) | (nbSecondNegativeReads>0))==TRUE)
      
      firstWin <- data.frame("Type"=rep("First",length(presentFirstWin)),
                             "Start" = presentFirstWin, 
                             "NbPositiveReads" = nbFirstPositiveReads[presentFirstWin], 
                             "NbNegativeReads" = nbFirstNegativeReads[presentFirstWin])
      secondWin <- data.frame("Type"=rep("Second",length(presentSecondWin)),
                              "Start" = presentSecondWin, "
                              NbPositiveReads" = nbSecondPositiveReads[presentSecondWin], 
                              "NbNegativeReads" = nbSecondNegativeReads[presentSecondWin])
      if (coverage){
        firstWin <- dplyr::mutate(firstWin,"MaxCoverage" = maxFirstCoverage)
        secondWin <- dplyr::mutate(secondWin,"MaxCoverage" = maxSecondCoverage)
      }
      ChromosomeFirst <- rep("",nrow(firstWin))
      ChromosomeSecond <- rep("",nrow(secondWin))
      for (i in seq_along(part)){
        id <- which(chromosomes == part[i])
        mi <- ceiling(statInfo$FirstBaseInPartition[id]/step)
        ma <- ceiling((statInfo$LastBaseInPartition[id]-win+1)/step)
        
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
  return(allWin %>% dplyr::mutate("Start" = (Start-1)*step+1))
}
