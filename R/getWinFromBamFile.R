#' @title get the strand information of all windows from bam files
#' @description get the number of positive/negative reads of all windows from bam files
#' @param files the input bam files. Your bamfiles should be sorted and have their index files located at the same path.
#' @param chromosomes the list of chromosomes to be read
#' @param mapqFilter every read that has mapping quality below \code{mapqFilter} will be removed before any analysis
#' @param partitionSize by default is 1e8, i.e. the bam file is read by blocks of chromosomes such that the total length of each block is at least 1e8
#' @param winWidth the width of the sliding window, 1000 by default.
#' @param winStep the step length to sliding the window, 100 by default.
#' @param paired if TRUE then the input bamfile will be considered as paired end reads. If missing, 100 thousands first reads will be inspected to test if the input bam file in paired end or single end.
#' @seealso filterDNA, summarizeHist, plotHist, plotWin
#' @export
#' @importFrom IRanges Views
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom Rsamtools bamMapqFilter<-
#' @importFrom Rsamtools bamMapqFilter
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools BamFileList
#' @examples
#' \dontrun{
#' file <- system.file("extdata","s1.sorted.bam",package = "strandCheckR")
#' win <- getWinFromBamFile(file)
#' }
getWinFromBamFile <- function(files, chromosomes, mapqFilter=0, partitionSize=1e8, winWidth=1000, winStep=100, paired){
  # Check the input is a BamFileList. Convert if necessary
  if (class(files) != "BamFileList") tryCatch(files <- BamFileList(files))
  
  # Check valid mapqFilter value
  if(mapqFilter < 0 || !is.numeric(mapqFilter)) stop("Invalid value for mapqFilter. Must be positive & numeric.")
  
  # Create a list to store the windows for each file
  allWin <- list(length(files))
  
  # Read through each bam file
  for (b in seq_along(files)){
    file <- files[[b]]
    
    # Check the paired status
    if (missing(paired)){ 
      paired <- checkPairedEnd(file$path)
    }
    
    # Get the seqinfo object & all genomic information
    sq <- seqinfo(file)
    allChromosomes <- seqnames(sq)
    
    # Subset chromosomes if required
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
    
    # Initialise the data frames for window information to be returned
    allWin[[b]] <- vector("list",length(partition))
    
    # Initialise the data frames containing the chromosome information in each partition
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
    if (mapqFilter) scanWhat <- c(scanWhat, "mapq")
    if (paired) scanWhat <- c(scanWhat,"flag")
    
    message("Reading file ",file$path)
    
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
      readInfo <- scanBam(file, param = sbp)
      nReadsInPart <- vapply(seq_along(readInfo),
                             function(i){length(readInfo[[i]]$strand)}, 
                             integer(1))
      chromosomeInfo$NbOriginalReads[idPart] <- nReadsInPart
      
      if (sum(nReadsInPart) > 0){
        
        # Calculate the first/last bases/reads in the chromosome partition
        chromosomeInfo[idPart,] <- chromosomeInfoInPartition(chromosomeInfo[idPart,], winStep)
        # concatenate several lists of the chromosome partition into one list
        readInfo <- concatenateAlignments(readInfo, chromosomeInfo[idPart,])
        if (paired){
          firstReadIndex <- ((floor(readInfo$flag/64) %% 2) == 1)
          secondReadIndex <- !firstReadIndex
          if (sum(firstReadIndex)==0){
            message("Only R2 reads are found")
            subset <- list(NULL) 
            type <- "R2"
          } else if (sum(secondReadIndex)==0){
            message("Only R1 reads are found")
            subset <- list(NULL)
            type <- "R1"
          } else {
            subset <- list("R1"=firstReadIndex,"R2"=secondReadIndex)
          }
        } else{
          subset <- list(NULL)
        }
        for (s in seq_along(subset)){
          
          winPositiveAlignments <- getWinOfAlignments(readInfo,"+",winWidth,winStep,readProp = 0, useCoverage=TRUE,subset[[s]])
          winNegativeAlignments <- getWinOfAlignments(readInfo,"-",winWidth,winStep,readProp = 0, useCoverage=TRUE,subset[[s]])
          
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
            allWin[[b]][[n]] <- win
            if (paired && length(subset)==0){
              allWin[[b]][[n]]$Type <- type
            }
          } 
          else{
            allWin[[b]][[n]]$Type <- Rle("R1")
            win$Type <- Rle("R2")
            allWin[[b]][[n]] <- rbind(allWin[[b]][[n]],win)
          }
        }
      }
      allWin[[b]][[n]]$Start <- (allWin[[b]][[n]]$Start-1)*winStep+1
    }
    # rbind all Partitions
    allWin[[b]] <- do.call(rbind, allWin[[b]])
    allWin[[b]]$File <- file$path
    allWin[[b]]$End <- allWin[[b]]$Start + winWidth
  }
  return(do.call(rbind,allWin))
}
