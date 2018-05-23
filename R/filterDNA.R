#' @title Filter Double Strand Sequences from a Bam File
#' @description Filter putative double strand DNA from a strand specific RNA-seq
#'  using a window sliding across the genome.
#' @param file the input bam file to be filterd. Your bamfile should be sorted
#'  and have an index file located at the same path.
#' @param mapqFilter every read that has mapping quality below \code{mapqFilter}
#'  will be removed before any analysis
#' @param destination The file path where the filtered output will be written
#' @param statfile the file to write the summary of the results
#' @param chromosomes the list of chromosomes to be filtered
#' @param partitionSize by default is 1e8, i.e. the bam file is read by block of
#'  chromosomes such that the total length of each block is at least 1e8
#' @param mustKeepRanges a GRanges object defines the ranges such that every
#'   read mapping to those ranges is always kept regardless the strand proportion
#'   of the windows containing them.
#' @param getWin if TRUE, the function will not only filter the bam file but
#'  also return a data frame containing the information of all windows. This
#'   data frame can also be obtained using \code{\link{getWinFromBamFile}}.
#'    FALSE by default.
#' @param winWidth the length of the sliding window, 1000 by default.
#' @param winStep the step length to sliding the window, 100 by default.
#' @param threshold if a window has strand proportion greater than \code{threshold},
#'  then the reads in that window will be kept. 0.7 by default
#' @param pvalueThreshold the threshold for the p-value in the test of keeping
#'  windows. 0.05 by default
#' @param useCoverage if TRUE, then the strand information in each window
#'  corresponds to the sum of coverage coming from positive/negative reads;
#'   and not the number of positive/negative reads as default.
#' @param min if \code{useCoverage=FALSE}, every window that has less than
#'  \code{min} reads will be rejected regardless the strand proportion. 
#'        If \code{useCoverage=TRUE}, every window has max coverage least
#'         than \code{min} will be rejected. 0 by default
#' @param max if \code{useCoverage=FALSE}, every window that has more than
#'  \code{max} reads will be kept regardless the strand proportion. 
#'  If \code{useCoverage=TRUE}, every window with coverage more than 
#'  \code{max} will be kept. 
#'  If 0 then it doesn't have effect on selecting window. 0 by default.
#' @param errorRate the probability that an RNA read takes the false strand.
#'  0.01 by default.
#' @param readProp A read is considered to be included in a window if more than
#'  \code{readProp} of it is in the window. 
#' Specified as a proportion.
#' @param paired if TRUE then the input bamfile will be considered as paired end
#'  reads. If missing, 100 thousands first reads will be inspected to test if
#'   the input bam file in paired end or single end.
#'
#' @details filterDNA reads a bam file containing strand specific RNA reads, and
#'  filter putative double strand DNA.
#' Using a window sliding across the genome, we calculate the positive/negative
#'  proportion of reads in each window.
#' We then use logistic regression to estimate the strand proportion of reads in
#'  each window, and calculate the p-value when compare that to the given threshold.
#' Let \eqn{\pi} be the strand proportion of reads in a window.
#' 
#' Null hypothesis for positive window: \eqn{\pi \le threshold}.
#' 
#' Null hypothesis for negative window: \eqn{\pi \ge 1-threshold}.
#' 
#' Only windows with p-value <= \code{pvalueThreshold} are kept. For a kept
#'  positive window, let P be its number of positive reads, and let M
#' be its number of negative reads. These M negative reads are supposed to come
#'  from double-strand DNA, then there should be also M postive reads among the
#' P positive reads come from double-strand DNA. In other words, there are only
#'  (P-M) positive reads come from RNA. Hence, each positive read in this window
#'   is kept
#' with the probability (P-M)/P. Each negative read is kept with the probability
#'  equalling the rate that an RNA read of your sample has wrong strand, which 
#'  is \code{errorRate}.
#' Since each alignment can be belonged to several windows, then the probability
#'  of keeping an alignment is the maximum probability defined by all windows
#'   that contain it.
#'
#' @seealso getWinFromBamFile, plotHist, plotWin
#'
#' @examples
#' \dontrun{
#' file <- system.file("extdata","s1.sorted.bam",package = "strandCheckR")
#' filterDNA(file,destination="out.bam",statfile = "out.stat")
#' }
#' 
#' @importFrom rbamtools bamReader getHeader bamWriter bamClose bamSave bamRange
#' @importFrom Rsamtools BamFile scanBam ScanBamParam
#' @importFrom GenomicRanges GRanges ranges
#' @importFrom GenomeInfoDb seqinfo seqnames seqlengths
#' @importFrom IRanges IRanges
#' @import S4Vectors
#'
#' @export
filterDNA <- function(file, destination, statfile, chromosomes, mapqFilter=0, 
                      partitionSize = 1e8, mustKeepRanges, getWin=FALSE, 
                      winWidth=1000, winStep=100, threshold=0.7, 
                      pvalueThreshold=0.05, min=0, max=0, errorRate=0.01, 
                      readProp=0.5, useCoverage=FALSE, paired){
  startTime <- proc.time()
  
  # Check the input is a BamFile. Convert if necessary
  if (length(file) > 1) {
    message("Multiple files provided. Only the first will be filtered")
    file <- file[1]
  }
  if (class(file) != "BamFile") tryCatch(file <- BamFile(file))
  
  # Check the destination path exists & has the suffix bam
  stopifnot(file.exists(dirname(destination)))
  destination <- ifelse(!grepl("\\.bam$", destination[1]), 
                        paste0(destination[1], ".bam"), 
                        destination[1])
  
  # Check valid mapqFilter value
  if(mapqFilter < 0 || !is.numeric(mapqFilter)) {
    stop("Invalid value for mapqFilter. Must be positive & numeric.")
  }
  
  # More input checks
  stopifnot(is.logical(getWin))
  winWidth <- tryCatch(as.integer(winWidth))
  winStep <- tryCatch(as.integer(winStep))
  stopifnot(threshold > 0 & threshold < 1)
  stopifnot(pvalueThreshold > 0 & pvalueThreshold < 1)
  stopifnot(all(is.numeric(max), is.numeric(min), is.numeric(errorRate), 
                is.numeric(readProp)))
  stopifnot(is.logical(useCoverage))
  
  # Check the paired status
  if (missing(paired)) paired <- checkPairedEnd(file$path)
  
  # Get the seqinfo object & all genomic information
  sq <- seqinfo(file)
  allRefChromosomes <- seqnames(sq)
  
  # Subset chromosomes if required
  if (!missing(chromosomes)){
    stopifnot(all(chromosomes %in% allRefChromosomes))
    sq <- sq[chromosomes]
  }
  else{
    chromosomes <- allRefChromosomes
  }
  nChr <- length(chromosomes)
  lengthSeq <- seqlengths(sq)
  
  # Allocate chromosomes for optimal partition sizes & speed
  partition <- partitionSeqinfo(sq, partitionSize = partitionSize)
  # Open a reader of the input bamfile to extract read afterward
  reader <- bamReader(file$path,idx=TRUE) 
  
  # Get the header of the input bam file and prepare a writer with the 
  # same header for the output bamfile.
  header <- getHeader(reader) 
  writer <- bamWriter(header,destination) 
  
  # Logit of the given threshold
  logitThreshold <- log(threshold/(1-threshold))
  
  # Define the file to write the summary of the results
  if (missing(statfile)){
    statfile <- "out.stat"
  }
  else{
    stopifnot(file.exists(dirname(statfile)))
  }
  message(paste0("Summary will be written to ",statfile))
  file.create(statfile)

  # Get the chromosome list of the ranges which must be kept
  if (!missing(mustKeepRanges)){
    stopifnot (class(mustKeepRanges) == "GRanges")
    allChromosomesMustKeep <- levels(seqnames(mustKeepRanges))
  }
  
  # Initialise the data frames containing sliding window information to be 
  # returned when getWin=TRUE
  if (getWin) allWin <- vector("list",length(partition))
  
  # Initialise the data frame containing the chromosome information in each partition
  chromosomeInfo <- data.frame(Sequence = chromosomes, 
                               Length = lengthSeq,
                               NbOriginalReads = rep(0, nChr),
                               NbKeptReads = rep(0, nChr),
                               FirstBaseInPartition = rep(NA, nChr),
                               LastBaseInPartition = rep(NA, nChr),
                               FirstReadInPartition = rep(NA, nChr),
                               LastReadInPartition = rep(NA, nChr),
                               stringsAsFactors = FALSE)
  
  # Define what to scan from bam file
  scanWhat <- c("pos", "cigar", "strand")
  if (mapqFilter) scanWhat <- c(scanWhat, "mapq")
  if (paired) scanWhat <- c(scanWhat, "flag")
  
  # Step through each part of partitions
  for (n in seq_along(partition)){
    
    # Set the required values for each partition
    part <- partition[[n]]
    idPart <- which(chromosomes %in% part)
    chrEnd <- lengthSeq[chromosomes %in% part]
    
    # Get & Summarise the reads
    # Return the reads from the bam file as a list, 
    # with each element containing reads from a single chr
    sbp <- ScanBamParam(what = scanWhat,
                        which = GRanges(seqnames = part,
                                        ranges = IRanges(start = 1,end = chrEnd)))
    readInfo <- scanBam(file, param = sbp)
    if (mapqFilter > 0){
      # Get only reads that pass the mapping quality filter
      # These need to be identified separately for the writing step
      mqfilter <- lapply(readInfo,function(chr){which(chr$mapq>=mapqFilter)})
      readInfo <- lapply(seq_along(readInfo),
                         function(i){chr <- lapply(scanWhat[scanWhat!="mapq"], 
                                                  function(field){field=readInfo[[i]][[field]][mqfilter[[i]]]})
                                    names(chr) <- scanWhat[scanWhat!="mapq"]
                                    chr})
      names(readInfo) <- names(mqfilter)
    } 
    nReadsInPart <- vapply(seq_along(readInfo),
                           function(i){length(readInfo[[i]]$strand)}, 
                           integer(1))
    chromosomeInfo$NbOriginalReads[idPart] <- nReadsInPart
    
    # Calculate the first/last bases/reads of each chromosome in the partition
    chromosomeInfo[idPart,] <- chromosomeInfoInPartition(chromosomeInfo[idPart,], winStep)
    
    # Initialise parameters
    keptAlignments <- c()

    if (sum(nReadsInPart) > 0){
      # Calculate the windows that overlap mustKeepRanges
      mustKeepWin <- list()
      if (!missing(mustKeepRanges)){
        if (length(intersect(allChromosomesMustKeep,part))>0){
          mustKeepWin <- getWinFromGranges(
            mustKeepRanges[seqnames(mustKeepRanges) %in% part], 
            chromosomeInfo[idPart,], winWidth, winStep)
        }
      }
      
      
      # Concatenate several lists of the chromosome partition into one list
      readInfo <- concatenateAlignments(readInfo,chromosomeInfo[idPart,])
      
      # Get the index of R1 and R2 reads and process each subset separately
      if (paired){
        firstReadIndex <- ((floor(readInfo$flag/64) %% 2) == 1)
        secondReadIndex <- !firstReadIndex
        if (sum(firstReadIndex)==0 || sum(secondReadIndex)==0) subset <- list(NULL) 
        else subset <- list("R1"=firstReadIndex,"R2"=secondReadIndex)
      } 
      else{
        subset <- list(NULL)
      }
      for (s in seq_along(subset)){
        
        # Get the ids of sliding windows containing each "+"/"-" read fragment
        winPositiveAlignments <- getWinOfAlignments(readInfo, "+", winWidth, 
                                                    winStep, readProp = readProp, 
                                                    useCoverage = (useCoverage || getWin), 
                                                    subset[[s]])
        winNegativeAlignments <- getWinOfAlignments(readInfo, "-", winWidth, 
                                                    winStep, readProp = readProp, 
                                                    useCoverage = (useCoverage || getWin), 
                                                    subset[[s]])
        
        # Calculate the keeping probability of each sliding window
        probaWin <- keptProbaWin(winPositiveAlignments, winNegativeAlignments, 
                                 winWidth, winStep, logitThreshold, pvalueThreshold,
                                 errorRate, mustKeepWin, min, max,getWin = getWin, 
                                 useCoverage=useCoverage) # the probability of each positive/negative window; this probability will be assigned to every positive/negative read in that window
        
        # Calculate the positive read fragments to be kept
        keptPositiveAlignment <- keptReadFragment(winPositiveAlignments$Win, 
                                                  probaWin$Positive, errorRate) 
        
        # Calculate the negative read fragments to be kept
        keptNegativeAlignment <- keptReadFragment(winNegativeAlignments$Win, 
                                                  probaWin$Negative, errorRate) 
       
        # If getWin=TRUE, then return the strand information of sliding windows
        if (getWin){
          win <- getWinInChromosome(probaWin$Win,part,chromosomeInfo[idPart,],
                                    winWidth, winStep)
          if (s==1){
            allWin[[n]] <- win
          } else{
            allWin[[n]]$Type <- Rle("R1")
            win$Type <- Rle("R2")
            allWin[[n]] <- rbind(allWin[[n]],win)
          }
        }

        # Infer the index of kept alignments within the partition
        kept <- c(unique(mcols(winPositiveAlignments$Win)$alignment[keptPositiveAlignment]),
                  unique(mcols(winNegativeAlignments$Win)$alignment[keptNegativeAlignment])) 
        # Add the current set to the vector of kept alignments
        keptAlignments <- sort(c(keptAlignments,kept))

        # Tidy up the memory a little
        rm(winPositiveAlignments,winNegativeAlignments)
        rm(probaWin)
      }
      rm(readInfo)
      
      # Write the kept alignments into the output bamfile
      for (i in seq_along(part)){#Run through each chromosome
        id <- which(chromosomes == part[i])
        
        if (chromosomeInfo$NbOriginalReads[id]>0){
          # Get the index within the partition of kept alignments coming from 
          # the current chromosome 
          index <- keptAlignments[(keptAlignments>=chromosomeInfo$FirstReadInPartition[id]) & 
                                    (keptAlignments<= chromosomeInfo$LastReadInPartition[id])]
           
          if (length(index)>0){
            chromosomeInfo$NbKeptReads[id] <- length(index)
            
            # Get the index within the chromosome of the kept alignments coming
            # from the current chromosome 
            index <- index - chromosomeInfo$FirstReadInPartition[id]+1
            
            # Get the index within the initial unfiltered chromosome of the kept
            # alignments coming from the current chromosome 
            if (mapqFilter>0) index <- mqfilter[[i]][index]
            
            # Write the kept alignments into the output bam file  
            bamSave(writer, bamRange(reader, 
                             c(which(allRefChromosomes == part[i]) - 1, 0, 
                               lengthSeq[id]))[index,])
          }
          rm(index)
        }
      }
      # Add the key information to the statfile
      cat(paste0("Sequence ",part,", length: ", chromosomeInfo$Length[idPart],
                 ", number of reads: ", chromosomeInfo$NbOriginalReads[idPart],
                 ", number of kept reads: ", chromosomeInfo$NbKeptReads[idPart],"\n"), 
          file=statfile, append = TRUE)
      rm(keptAlignments)
    }
  }
  bamClose(writer)
  bamClose(reader)
  
  cat("Summary:\n", file = statfile, append = TRUE)
  cat(paste0("Number of original reads: ", sum(chromosomeInfo$NbOriginalReads),
             ", number of kept reads: ", sum(chromosomeInfo$NbKeptReads),
             ", removal proportion: ",
             (sum(chromosomeInfo$NbOriginalReads)-sum(chromosomeInfo$NbKeptReads))/sum(chromosomeInfo$NbOriginalReads)), 
      "\n", file = statfile, append = TRUE)
  endTime <- proc.time()
  cat(paste0("Total elapsed time: ", (endTime-startTime)[[3]]/60," minutes\n"), 
      file = statfile, append = TRUE)
  if (getWin){
    return(do.call(rbind,allWin))
  }
}



