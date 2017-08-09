#' @title get the number of positive/negative reads of all windows from a single end bam file
#'
#' @param file the input single-end bam file. Your bamfile should be sorted and have an index file located at the same path as well.
#' @param chromosomes the list of chromosomes to be read
#' @param mapqFilter very read that has mapping quality below \code{mapqFilter} will be removed before any analysis
#' @param partitionSize by default is 1e8, i.e. the bam file is read by blocks of chromosomes 
#' such that the total length of each block is at least 1e8
#' @param winWidth the width of the sliding window, 1000 by default.
#' @param winStep the winStep length to sliding the window, 100 by default.
#' 
#' @seealso filterDNA, filterDNAPairs, getWinFromPairedBamFile, plotHist, plotWin
#' @export
#' @importFrom IRanges Views
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom Rsamtools bamMapqFilter
#' @importFrom Rsamtools ScanBamParam
#' @examples
#' file <- system.file("data","s1.chr1.bam",package = "strandCheckR")
#' win <- getWinFromBamFile(file)
#'
getWinFromBamFile <- function(file, chromosomes, mapqFilter=0, partitionSize=1e8, winWidth=1000, winStep=100, paired){
  
  # Check the input is a BamFile. Convert if necessary
  if (class(file) != "BamFile") file <- BamFile(file)
  
  # Check valid mapqFilter value
  if(mapqFilter < 0 || !is.numeric(mapqFilter)) stop("Invalid value for mapqFilter. Must be positive & numeric.")
  
  # Check the paired status
  if (missing(paired)){
    #autodetect
  }
  
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
  allWin <- list(length(partition))
  statInfo <- data.frame(Sequence = chromosomes, 
                         Length = lengthSeq,
                         NbOriginalReads = rep(NA, nChr),
                         FirstBaseInPartition = rep(NA, nChr),
                         LastBaseInPartition = rep(NA, nChr),
                         FirstReadInPartition = rep(NA,nChr),
                         LastReadInPartition = rep(NA, nChr),
                         stringsAsFactors = FALSE)
  
  # Step through each set of partitions
  for (n in seq_along(partition)){
    
    # Set the required values for this section
    part <- partition[[n]]
    idPart <- which(chromosomes %in% part)
    chrEnd <- lengthSeq[chromosomes %in% part]
    
    # Get & Summarise the reads
    sbp <- ScanBamParam(what = c("pos","cigar","strand"),
                        which = GRanges(seqnames = part,ranges = IRanges(start = 1,end = chrEnd)))
    bamMapqFilter(sbp) <- mapqFilter
    
    # Return the reads from the bam file as a list, with each element containing reads from a single chr
    bam <- scanBam(file, param = sbp)
    nReadsInPart <- vapply(seq_along(bam),
                           function(i){length(bam[[i]]$strand)}, 
                           integer(1))
    statInfo$NbOriginalReads[idPart] <- nReadsInPart
    
    if (sum(nReadsInPart) > 0){
      
      # Calculate the first/last bases/reads in the chromosome partition
      statInfo[idPart,] <- statInfoInPartition(statInfo[idPart,], winStep)
      # concatenate several lists of the chromosome partition into one list
      bam <- concatenateAlignments(bam, statInfo[idPart,])
      winPositiveAlignments <- getWinOfAlignments(bam, "+", winWidth, winStep, limit = 1, coverage=TRUE)
      winNegativeAlignments <- getWinOfAlignments(bam, "-", winWidth, winStep, limit = 1, coverage=TRUE)
      rm(bam)
      
      ##################################################
      # calculate strand information based on coverage #
      ##################################################
      
      # make sure winPositiveAlignments$Coverage and winNegativeAlignments$Coverage
      # have the same length to avoid some warnings afterward
      lastBase <- max(c(length(winPositiveAlignments$Coverage),
                        length(winNegativeAlignments$Coverage)))
      lenPC <- length(winPositiveAlignments$Coverage)
      lenNC <- length(winNegativeAlignments$Coverage)
      if (lenNC < lastBase) {
        winNegativeAlignments$Coverage <- c(winNegativeAlignments$Coverage,rep(0, lastBase - lenNC))
      } 
      if (lenPC < lastBase){
        winPositiveAlignments$Coverage <- c(winPositiveAlignments$Coverage,rep(0, lastBase - lenPC))
      }
      
      #calculate the number of positive and negative bases in each window
      nbWin <- ceiling((lastBase - winWidth) / winStep) + 1
      st <- seq(1, (nbWin - 1)*winStep + 1, by =winStep)
      end <- seq(winWidth, (nbWin-1)*winStep + winWidth, by = winStep)
      CovPositive <- Views(winPositiveAlignments$Coverage, start = st, end = end) 
      CovPositive <- Rle(sum(CovPositive))
      CovNegative <- Views(winNegativeAlignments$Coverage, start = st, end = end) 
      CovNegative <- Rle(sum(CovNegative))
      
      #calculate max the max coverage in each window
      maxCoverage <- Rle(max(Views(winPositiveAlignments$Coverage + winNegativeAlignments$Coverage, 
                                   start = st, end=end)))
      
      ######################################################
      # calculate strand information based on nbr of reads #
      ######################################################
      
      # Calculate strand information based on number of reads 
      # have the same length to avoid some warnings afterward
      NbPositive <- coverage(winPositiveAlignments$Win)
      NbNegative <- coverage(winNegativeAlignments$Win)
      # Find the last window in both sets of windows
      lastWin <- max(c(end(winNegativeAlignments$Win), 
                       end(winPositiveAlignments$Win)))
      stopifnot(lastWin == nbWin)
      # Fill with zeroes if required
      #make sure NbPositive and NbNegative have the same length
      lenP <- length(NbPositive)
      lenN <- length(NbNegative)
      if (lenN < lastWin) NbNegative <- c(NbNegative,rep(0,lastWin - lenN))
      if (lenP < lastWin) NbPositive <- c(NbPositive,rep(0,lastWin - lenP))
      
      
      #fill the information of the present window into the data frame to be returned
      presentWin <- which(as.vector( (NbPositive>0) | (NbNegative>0))==TRUE)
      Chromosome <- Rle(rep("", length(presentWin)))
      Start <- integer(length(presentWin))
      # Step through for each chromosome in the current partition
      for (i in seq_along(part)){
        currentChr <- part[i]
        id <- which(chromosomes == currentChr)
        idFirst <- ceiling(statInfo$FirstBaseInPartition[id] / winStep) # id of the first window of the chromosome
        idLast <- ceiling((statInfo$LastBaseInPartition[id] - winWidth+1) / winStep)# id of the last window of the chromosome
        idRows <- which(presentWin >= idFirst & presentWin <= idLast) #get the windows of the chromosome
        Chromosome[idRows] <- currentChr
        Start[idRows] <- (presentWin[idRows] - idFirst)*winStep + 1
      }
      
      Win <- DataFrame(Chr = Chromosome, Start = Start, 
                       NbPositive = NbPositive[presentWin], NbNegative = NbNegative[presentWin],
                       CovPositive = CovPositive[presentWin], CovNegative = CovNegative[presentWin],
                       MaxCoverage = maxCoverage[presentWin])
      
      allWin[[n]] <- Win
      
    }
  }
  
  # rbind all Partitions
  do.call(rbind, allWin)
  
}
