#' @title get the strand information of all windows from a single end bam file
#' @description get the number of positive/negative reads of all windows from a single end bam file
#' @param file the input single-end bam file. Your bamfile should be sorted and have an index file located at the same path as well.
#' @param chromosomes the list of chromosomes to be read
#' @param mapqFilter very read that has mapping quality below \code{mapqFilter} will be removed before any analysis
#' @param partitionSize by default is 1e8, i.e. the bam file is read by blocks of chromosomes 
#' such that the total length of each block is at least 1e8
#' @param winWidth the width of the sliding window, 1000 by default.
#' @param winStep the winStep length to sliding the window, 100 by default.
#' @param paired if TRUE then the input bamfile will be considered as paired end reads. If missing, 100 thousands first reads will be inspected to test if the input bam file in paired end or single end.
#' @seealso filterDNA, getWinFromPairedBamFile, plotHist, plotWin
#' @export
#' @importFrom IRanges Views
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom Rsamtools bamMapqFilter<-
#' @importFrom Rsamtools ScanBamParam
#' @examples
#' file <- system.file("extdata","s1.chr1.bam",package = "strandCheckR")
#' win <- getWinFromBamFile(file)
#'
getWinFromBamFile <- function(file, chromosomes, mapqFilter=0, partitionSize=1e8, winWidth=1000, winStep=100, paired){
  # Check the input is a BamFile. Convert if necessary
  if (class(file) != "BamFile") file <- BamFile(file)
  
  # Check valid mapqFilter value
  if(mapqFilter < 0 || !is.numeric(mapqFilter)) stop("Invalid value for mapqFilter. Must be positive & numeric.")
  
  # Check the paired status
  if (missing(paired)){
    message("Testing paired end bam file by checking the first 100000 reads")
    checkFile <- BamFile(file$path, yieldSize = 100000)
    flag <- scanBam(checkFile, param = ScanBamParam(what = "flag"))[[1]]$flag
    paired <- any(flag %% 2 == 1)
    if (paired) message("Your bam file is paired end") else message("Your bam file is singe end")
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
  allWin <- vector("list",length(partition))
  chromosomeInfo <- data.frame(Sequence = chromosomes, 
                               Length = lengthSeq,
                               NbOriginalReads = rep(NA, nChr),
                               FirstBaseInPartition = rep(NA, nChr),
                               LastBaseInPartition = rep(NA, nChr),
                               FirstReadInPartition = rep(NA, nChr),
                               LastReadInPartition = rep(NA, nChr),
                               stringsAsFactors = FALSE)
  
  #what to scan from bam file
  scanWhat <- c("pos","cigar","strand")
  if (paired) scanWhat <- c(scanWhat,"flag")
  # Step through each set of partitions
  for (n in seq_along(partition)){
    
    # Set the required values for this section
    part <- partition[[n]]
    idPart <- which(chromosomes %in% part)
    chrEnd <- lengthSeq[chromosomes %in% part]
    
    # Get & Summarise the reads
    sbp <- ScanBamParam(what = scanWhat,
                          which = GRanges(seqnames = part,ranges = IRanges(start = 1,end = chrEnd)))
    bamMapqFilter(sbp) <- mapqFilter
    
    # Return the reads from the bam file as a list, with each element containing reads from a single chr
    bam <- scanBam(file, param = sbp)
    nReadsInPart <- vapply(seq_along(bam),
                           function(i){length(bam[[i]]$strand)}, 
                           integer(1))
    chromosomeInfo$NbOriginalReads[idPart] <- nReadsInPart
    
    if (sum(nReadsInPart) > 0){
      
      # Calculate the first/last bases/reads in the chromosome partition
      chromosomeInfo[idPart,] <- chromosomeInfoInPartition(chromosomeInfo[idPart,], winStep)
      # concatenate several lists of the chromosome partition into one list
      bam <- concatenateAlignments(bam, chromosomeInfo[idPart,])
      if (paired){
        firstReadIndex <- ((floor(bam$flag/64) %% 2) == 1)
        secondReadIndex <- !firstReadIndex
        subset <- list("first"=firstReadIndex,"second"=secondReadIndex)
      } else{
        subset <- list(NULL)
      }
      for (s in seq_along(subset)){
        winPositiveAlignments <- getWinOfAlignments(bam,"+",winWidth,winStep,limit = 0, useCoverage=TRUE,subset[[s]])
        winNegativeAlignments <- getWinOfAlignments(bam,"-",winWidth,winStep,limit = 0, useCoverage=TRUE,subset[[s]])
        
        ##################################################
        # calculate strand information based on coverage #
        ##################################################
        fromCoverage <- calculateStrandCoverage(winPositiveAlignments,winNegativeAlignments,winWidth,winStep)
        
        ######################################################
        # calculate strand information based on nbr of reads #
        ######################################################
        fromNbReads <- calculateStrandNbReads(winPositiveAlignments,winNegativeAlignments)
        
        stopifnot(length(fromCoverage$CovPositive) == length(fromNbReads$NbPositive))
        
        #fill the information of the present window into the data frame to be returned
        presentWin <- which(as.vector( (fromNbReads$NbPositive>0) | (fromNbReads$NbNegative>0))==TRUE)
        win <- DataFrame(Chr = Rle(part[1]),
                         Start = presentWin, 
                         NbPositive = fromNbReads$NbPositive[presentWin], NbNegative = fromNbReads$NbNegative[presentWin],
                         CovPositive = fromCoverage$CovPositive[presentWin], CovNegative = fromCoverage$CovNegative[presentWin],
                         MaxCoverage = fromCoverage$MaxCoverage[presentWin])
        win <- getWinInChromosome(win,part,chromosomeInfo[idPart,],winWidth,winStep)
        if (s==1){
          allWin[[n]] <- win
        } else{
          allWin[[n]]$Type <- Rle("First")
          win$Type <- Rle("Second")
          allWin[[n]] <- rbind(allWin[[n]],win)
        }
      }
    }
  }
  # rbind all Partitions
  do.call(rbind, allWin)
}
