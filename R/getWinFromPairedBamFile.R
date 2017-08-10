#' @title get the number of positive/negatives of all windows from a paired end bam file
#'
#' @param bamfilein the input paired-end bam file. Your bamfile should be sorted and have an index file located at the same path as well.
#' @param mq every read that has mapping quality below \code{mq} will be removed before any analysis
#' @param chromosomes the list of chromosomes to be read
#' @param partitionSize by default is 1e8, i.e. the bam file is read by block of chromosomes such that the total length of each block is at least 1e8
#' @param winWidth the length of the sliding window, 1000 by default.
#' @param winStep the step length to sliding the window, 100 by default.
#' @param limit a read is considered to be included in a window if and only if at least \code{limit} percent of it is in the window. 0.75 by default
#' @param coverage if TRUE, then the strand information in each window corresponds to the sum of coverage coming from positive/negative reads; and not the number of positive/negative reads as default.
#'
#' @seealso filterDNA, filterDNAPairs, getWinFromBamFile, plotHist, plotWin
#' @export
#' @importFrom IRanges Views
#' @importFrom dplyr mutate
#' @examples#' 
#' bamfilein <- system.file("data","120.10.bam",package = "rnaCleanR")
#' win <- getWinFromPairedBamFile(bamfilein)
#'
getWinFromPairedBamFile <- function(file, chromosomes, mapqFilter=0, partitionSize=1e8, winWidth=1000,winStep=100,limit=0.75,coverage=FALSE){
  # Check the input is a BamFile. Convert if necessary
  if (class(file) != "BamFile") file <- BamFile(file)
  
  # Check valid mapqFilter value
  if(mapqFilter < 0 || !is.numeric(mapqFilter)) stop("Invalid value for mapqFilter. Must be positive & numeric.")
  
  # Get the seqinfo object & all genomic information
  sq <- seqinfo(file)
  allChromosomes <- seqnames(sq)
  # Subset if required
  if (!missing(chromosomes)){
    stopifnot(all(chromosomes %in% allChromosomes))
    sq <- sq[chromosomes]
  }
  else{
    chromosomes <- allChromosomes
  }
  nChr <- length(chromosomes)
  lengthSeq <- seqlengths(sq)
  
  
  # Allocate chromosomes for optimal partition sizes & speed
  partition <- partitionSeqinfo(sq, partitionSize = partitionSize)
  
  # Initialise the data frames for window information 
  allWin <- vector("list",length(partition))
  
  statInfo <- data.frame(Sequence = chromosomes, 
                         Length = lengthSeq,
                         "NbOriginalReads" = rep(NA,length(chromosomes)), 
                         "NbOriginalFirstReads" = rep(NA,length(chromosomes)), 
                         "NbOriginalSecondReads" = rep(NA,length(chromosomes)), 
                         "FirstBaseInPartition" = rep(NA,length(chromosomes)),
                         "LastBaseInPartition" = rep(NA,length(chromosomes)),
                         "FirstReadInPartition" = rep(NA,length(chromosomes)),
                         "LastReadInPartition" = rep(NA,length(chromosomes)),
                         stringsAsFactors = FALSE)
  
  for (n in seq_along(partition)){
    # Set the required values for this section
    part <- partition[[n]]
    idPart <- which(chromosomes %in% part)
    chrEnd <- lengthSeq[chromosomes %in% part]
    
    # Get & Summarise the reads
    sbp <- ScanBamParam(what = c("pos","cigar","strand","flag"),
                        which = GRanges(seqnames = part,ranges = IRanges(start = 1,end = chrEnd)))
    bamMapqFilter(sbp) <- mapqFilter
    
    # Return the reads from the bam file as a list, with each element containing reads from a single chr
    bam <- scanBam(file, param = sbp)
    nReadsInPart <- vapply(seq_along(bam),
                           function(i){length(bam[[i]]$strand)}, 
                           integer(1))
    statInfo$NbOriginalReads[idPart] <- nReadsInPart
    
    if (sum(statInfo$NbOriginalReads[idPart])>0){
      statInfo[idPart,] <- statInfoInPartition(statInfo[idPart,],winStep)
      bam <- concatenateAlignments(bam, statInfo[idPart,])
      firstReadIndex <- ((floor(bam$flag/64) %% 2) == 1)
      secondReadIndex <- !firstReadIndex
      winFirstPositiveAlignments <- getWinOfAlignments(bam,"+",winWidth,winStep,limit = 1, useCoverage=TRUE,firstReadIndex)
      winSecondPositiveAlignments <- getWinOfAlignments(bam,"+",winWidth,winStep,limit = 1, useCoverage=TRUE,secondReadIndex)
      winFirstNegativeAlignments <- getWinOfAlignments(bam,"-",winWidth,winStep,limit = 1, useCoverage=TRUE,firstReadIndex)
      winSecondNegativeAlignments <- getWinOfAlignments(bam,"-",winWidth,winStep,limit = 1, useCoverage=TRUE,secondReadIndex)
      
      ##################################################
      # calculate strand information based on coverage #
      ##################################################
      
      lenPC <- length(winFirstPositiveAlignments$Coverage)
      lenNC <- length(winFirstNegativeAlignments$Coverage)
      lastBase <- max(lenPC,lenNC)
      
      if (lenNC < lastBase) {
        winFirstNegativeAlignments$Coverage <- c(winFirstNegativeAlignments$Coverage,rep(0,lastBase-lenNC))
      }
      if (lenPC < lastBase) {
        winFirstPositiveAlignments$Coverage <- c(winFirstPositiveAlignments$Coverage,rep(0,lastBase-lenPC))
      }
      nbWin <- ceiling((lastBase-winWidth)/winStep)+1
      CovFirstPositive <- Views(winFirstPositiveAlignments$Coverage,
                                start = seq(1,(nbWin-1)*winStep+1,winStep),
                                end=seq(winWidth,(nbWin-1)*winStep+winWidth,winStep)) 
      CovFirstPositive <- Rle(sum(CovFirstPositive))
      CovFirstNegative <- Views(winFirstNegativeAlignments$Coverage,
                                start = seq(1,(nbWin-1)*winStep+1,winStep),
                                end=seq(winWidth,(nbWin-1)*winStep+winWidth,winStep)) 
      CovFirstNegative <- Rle(sum(CovFirstNegative))
      MaxFirstCoverage <- Rle(max(Views(winFirstPositiveAlignments$Coverage+winFirstNegativeAlignments$Coverage,
                                        start = seq(1,(nbWin-1)*winStep+1,winStep),
                                        end=seq(winWidth,(nbWin-1)*winStep+winWidth,winStep))))
      
      lenPC <- length(winSecondPositiveAlignments$Coverage)
      lenNC <- length(winSecondNegativeAlignments$Coverage)
      lastBase <- max(lenPC,lenNC)
      nbWin <- ceiling((lastBase-winWidth)/winStep)+1
      if (lenNC < lastBase) {
        winSecondNegativeAlignments$Coverage <- c(winSecondNegativeAlignments$Coverage,rep(0,lastBase-lenNC))
      } 
      if (lenPC < lastBase) {
        winSecondPositiveAlignments$Coverage <- c(winSecondPositiveAlignments$Coverage,rep(0,lastBase-lenPC))
      }
      CovSecondPositive <- Views(winSecondPositiveAlignments$Coverage,
                                 start = seq(1,(nbWin-1)*winStep+1,winStep),
                                 end=seq(winWidth,(nbWin-1)*winStep+winWidth,winStep)) 
      CovSecondPositive <- Rle(sum(CovSecondPositive))
      CovSecondNegative <- Views(winSecondNegativeAlignments$Coverage,
                                 start = seq(1,(nbWin-1)*winStep+1,winStep),
                                 end=seq(winWidth,(nbWin-1)*winStep+winWidth,winStep))
      CovSecondNegative <- Rle(sum(CovSecondNegative))
      MaxSecondCoverage <- Rle(max(Views(winSecondPositiveAlignments$Coverage+winSecondNegativeAlignments$Coverage,
                                         start = seq(1,(nbWin-1)*winStep+1,winStep),
                                         end=seq(winWidth,(nbWin-1)*winStep+winWidth,winStep))))
      
      ######################################################
      # calculate strand information based on nbr of reads #
      ######################################################
      
      NbFirstPositive <- coverage(winFirstPositiveAlignments$Win)
      NbSecondPositive <- coverage(winSecondPositiveAlignments$Win)
      NbFirstNegative <- coverage(winFirstNegativeAlignments$Win)
      NbSecondNegative <- coverage(winSecondNegativeAlignments$Win)
      lenFirstP <- length(NbFirstPositive)
      lenSecondP <- length(NbSecondPositive)
      lenFirstN <- length(NbFirstNegative)
      lenSecondN <- length(NbSecondNegative)
      lastWinFirst <- max(lenFirstP,lenFirstN)
      lastWinSecond <- max(lenSecondP,lenSecondN)
      if (lenFirstN < lastWinFirst) {
        NbFirstNegative <- c(NbFirstNegative,rep(0,lastWinFirst - lenFirstN))
      } 
      if (lenFirstP < lastWinFirst) {
        NbFirstPositive <- c(NbFirstPositive,rep(0,lastWinFirst - lenFirstP))
      }
      if (lenSecondN < lastWinSecond) {
        NbSecondNegative <- c(NbSecondNegative,rep(0,lastWinSecond - lenSecondN))
      } 
      if (lenSecondP < lastWinSecond) {
        NbSecondPositive <- c(NbSecondPositive,rep(0,lastWinSecond - lenSecondP))
      }
      presentFirstWin <- which(as.vector((NbFirstPositive>0) | (NbFirstNegative>0))==TRUE)
      presentSecondWin <- which(as.vector((NbSecondPositive>0) | (NbSecondNegative>0))==TRUE)
      
      ChromosomeFirst <- Rle(rep("", length(presentFirstWin)))
      ChromosomeSecond <- Rle(rep("", length(presentSecondWin)))
      StartFirst <- integer(length(presentFirstWin))
      StartSecond <- integer(length(presentSecondWin))
      # Step through for each chromosome in the current partition
      for (i in seq_along(part)){
        currentChr <- part[i]
        id <- which(chromosomes == currentChr)
        idFirst <- ceiling(statInfo$FirstBaseInPartition[id] / winStep) # id of the first window of the chromosome
        idLast <- ceiling((statInfo$LastBaseInPartition[id] - winWidth+1) / winStep)# id of the last window of the chromosome
        idRowsFirst <- which(presentFirstWin >= idFirst & presentFirstWin <= idLast) #get the windows of the chromosome
        idRowsSecond <- which(presentSecondWin >= idFirst & presentSecondWin <= idLast) #get the windows of the chromosome
        ChromosomeFirst[idRowsFirst] <- currentChr
        ChromosomeSecond[idRowsSecond] <- currentChr
        StartFirst[idRowsFirst] <- (presentFirstWin[idRowsFirst] - idFirst)*winStep + 1
        StartSecond[idRowsSecond] <- (presentSecondWin[idRowsSecond] - idFirst)*winStep + 1
      }
      allWin[[n]] <- rbind(DataFrame(Type="First",Chr = ChromosomeFirst, Start = StartFirst, 
                                     NbPositive = NbFirstPositive[presentFirstWin], NbNegative = NbFirstNegative[presentFirstWin],
                                     CovPositive = CovFirstPositive[presentFirstWin], CovNegative = CovFirstNegative[presentFirstWin],
                                     MaxCoverage = MaxFirstCoverage[presentFirstWin]),
                           DataFrame(Type="Second",Chr = ChromosomeSecond, Start = StartSecond, 
                                     NbPositive = NbSecondPositive[presentSecondWin], NbNegative = NbSecondNegative[presentSecondWin],
                                     CovPositive = CovSecondPositive[presentSecondWin], CovNegative = CovSecondNegative[presentSecondWin],
                                     MaxCoverage = MaxSecondCoverage[presentSecondWin]))
    }
  }
  do.call(rbind,allWin)
}
