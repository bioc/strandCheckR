#' @title get the strand information of all windows from a paired end bam file
#' @description get the number of positive/negative reads of all windows from a paired end bam file
#' @param file the input paired-end bam file. Your bamfile should be sorted and have an index file located at the same path as well.
#' @param mapqFilter every read that has mapping quality below \code{mapqFilter} will be removed before any analysis
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
#' bamfilein <- system.file("data","120.10.bam",package = "strandCheckR")
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
      
      fromCoverageFirst <- calculateStrandCoverage(winFirstPositiveAlignments,winFirstNegativeAlignments,winWidth,winStep)
      fromCoverageSecond <- calculateStrandCoverage(winSecondPositiveAlignments,winSecondNegativeAlignments,winWidth,winStep)
      
      
      ######################################################
      # calculate strand information based on nbr of reads #
      ######################################################
      
      fromNbReadsFirst <- calculateStrandNbReads(winFirstPositiveAlignments,winFirstNegativeAlignments)
      fromNbReadsSecond <- calculateStrandNbReads(winSecondPositiveAlignments,winSecondNegativeAlignments)
      
      
      stopifnot(length(fromCoverageFirst$CovPositive) == length(fromNbReadsFirst$NbPositive))
      stopifnot(length(fromCoverageSecond$CovPositive) == length(fromNbReadsSecond$NbPositive))
      
      #fill the information of the present window into the data frame to be returned
      presentWinFirst <- which(as.vector( (fromNbReadsFirst$NbPositive>0) | (fromNbReadsFirst$NbNegative>0))==TRUE)
      firstWin <- DataFrame(Start = presentWinFirst, 
                               NbPositive = fromNbReadsFirst$NbPositive[presentWinFirst], NbNegative = fromNbReadsFirst$NbNegative[presentWinFirst],
                               CovPositive = fromCoverageFirst$CovPositive[presentWinFirst], CovNegative = fromCoverageFirst$CovNegative[presentWinFirst],
                               MaxCoverage = fromCoverageFirst$MaxCoverage[presentWinFirst])
      firstWin <- getWinInChromosome(firstWin,part,statInfo[idPart,],winWidth,winStep)
      firstWin$Type <- "First"
      
      presentWinSecond <- which(as.vector( (fromNbReadsSecond$NbPositive>0) | (fromNbReadsSecond$NbNegative>0))==TRUE)
      secondWin <- DataFrame(Start = presentWinSecond, 
                             NbPositive = fromNbReadsSecond$NbPositive[presentWinSecond], NbNegative = fromNbReadsSecond$NbNegative[presentWinSecond],
                             CovPositive = fromCoverageSecond$CovPositive[presentWinSecond], CovNegative = fromCoverageSecond$CovNegative[presentWinSecond],
                             MaxCoverage = fromCoverageSecond$MaxCoverage[presentWinSecond])
      secondWin <- getWinInChromosome(secondWin,part,statInfo[idPart,],winWidth,winStep)
      secondWin$Type <- "Second"
      allWin[[n]] <- rbind(firstWin,secondWin)
    }
  }
  do.call(rbind,allWin)
}
