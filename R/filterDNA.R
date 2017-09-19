#' @title Filter Single End Bam File
#'
#' @description Filter putative double strand DNA from a strand specific single-end RNA-seq using a window sliding across the genome.

#'
#' @param file the input single end bam file to be filterd. Your bamfile should be sorted and have an index file located at the same path as well.
#' @param mapqFilter every read that has mapping quality below \code{mapqFilter} will be removed before any analysis
#' @param fileout the output filtered bam file
#' @param statfile the file to write the summary of the results
#' @param chromosomes the list of chromosomes to be filtered
#' @param partitionSize by default is 1e8, i.e. the bam file is read by block of chromosomes such that the total length of each block is at least 1e8
#' @param mustKeepRanges a GRanges object defines the ranges such that every read maps to those ranges must be always kept regardless the strand proportion of the windows containing them.
#' @param getWin if TRUE, the function will return a data frame containing the information of all windows. It's FALSE by default.
#' @param winWidth the length of the sliding window, 1000 by default.
#' @param winStep the step length to sliding the window, 100 by default.
#' @param threshold the threshold upper which we keep the windows. 0.7 by default
#' @param pvalueThreshold the threshold for the p-value. 0.05 by default
#' @param min In the case that \code{useCoverage=FALSE}, if a window has least than \code{min} reads, then it will be rejected regardless the strand proportion. 
#'        For the case that \code{useCoverage=TRUE}, if a window has max coverage least than \code{min}, then it will be rejected. 0 by default
#' @param max In the case that \code{useCoverage=FALSE},if a window has more than \code{max} reads, then it will be kept regardless the strand proportion. 
#'        For the case that \code{useCoverage=TRUE}, if a window has max coverage more than \code{max}, then it will be kept. 
#'        If 0 then it doesn't have effect on selecting window. 0 by default.
#' @param errorRate the probability that an RNA read takes the false strand. 0.01 by default
#' @param readProp a read is considered to be included in a window if and only if at least \code{readProp} percent of it is in the window. 0.75 by default
#' @param useCoverage if TRUE, then the strand information in each window corresponds to the sum of coverage coming from positive/negative reads; and not the number of positive/negative reads as default.
#' @param paired if TRUE then the input bamfile will be considered as paired end reads. If missing, 100 thousands first reads will be inspected to test if the input bam file in paired end or single end.
#'
#' @details filterDNA reads a single end bam file containing strand specific RNA reads, and filter putative double strand DNA.
#' Using a window sliding across the genome, we calculate the positive/negative proportion of reads in that window.
#' For each window, we use logistic regression to estimate the proportion of reads in the window derived from
#' stranded RNA (positive or negative).
#'
#' Let \eqn{\pi} be the proportion of reads in the window derived from stranded RNA (positive or negative)
#'
#' Null hypothesis: \eqn{\pi < {\pi}_{0}} where \eqn{{\pi}_{0}} is the given threshold.
#'
#' Only windows with p-value <= 0.05 are kept. For a positive window that is kept, let P be its number of positive reads, and let M
#' be its number of negative reads. Since these M negative reads should come from double-strand DNA, then there should be also M postive reads among the
#' P positive reads come from double-strand DNA. In other words, there are only (P-M) positive reads come from RNA. Hence, each positive read in this window is kept
#' with the probability equalling (P-M)/P. Each negative read is kept with the probability equalling the given \code{errorRate} which is the rate that
#' an RNA read of your sample has wrong strand.
#'
#' Since each alignment can be belonged to several windows, then the probability of keeping an alignment is the maximum probability defined by
#' all windows that contain it.
#'
#' @seealso getWinFromBamFile, plotHist, plotWin
#'
#' @examples
#' \dontrun{
#' file <- system.file("extdata","s1.chr1.bam",package = "strandCheckR")
#' filterDNA(file,fileout="out.bam",statfile = "out.stat")
#' }
#' 
#' @importFrom rbamtools bamReader getHeader bamWriter bamClose bamSave bamRange
#' @importFrom Rsamtools BamFile scanBam ScanBamParam
#' @importFrom GenomicRanges GRanges ranges
#' @importFrom GenomeInfoDb seqinfo seqnames seqlengths
#' @importFrom IRanges IRanges
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @import S4Vectors
#'
#' @export
filterDNA <- function(file,fileout,statfile,chromosomes,mapqFilter=0,partitionSize = 1e8, mustKeepRanges,getWin=FALSE,winWidth=1000,winStep=100,threshold=0.7,pvalueThreshold=0.05,min=0,max=0,errorRate=0.01,readProp=0.75,useCoverage=FALSE,paired){
  startTime <- proc.time()
  
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
  allRefChromosomes <- seqnames(sq)
  # Subset if required
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
  
  reader <- bamReader(file$path,idx=TRUE) #open a reader of the input bamfile to extract read afterward
  header <- getHeader(reader) #get the header of the input bam file
  writer <- bamWriter(header,fileout) #prepare to write the output bamfile with the same header
  logitThreshold <- log(threshold/(1-threshold))
  
  if (missing(statfile)){
    message("Summary will be written to file out.stat")
    statfile <- "out.stat"
  }

  if (!missing(mustKeepRanges)){
    allChromosomesMustKeep <- levels(seqnames(mustKeepRanges))
  }
  
  # Initialise the data frames for window information 
  allWin <- vector("list",length(partition))
  chromosomeInfo <- data.frame(Sequence = chromosomes, 
                               Length = lengthSeq,
                               NbOriginalReads = rep(0, nChr),
                               NbKeptReads = rep(0, nChr),
                               FirstBaseInPartition = rep(NA, nChr),
                               LastBaseInPartition = rep(NA, nChr),
                               FirstReadInPartition = rep(NA, nChr),
                               LastReadInPartition = rep(NA, nChr),
                               stringsAsFactors = FALSE)
  append <- FALSE
  
  #what to scan from bam file
  scanWhat <- c("pos","cigar","strand")
  if (mapqFilter) scanWhat <- c(scanWhat,"mapq")
  if (paired) scanWhat <- c(scanWhat,"flag")
  
  # Step through each set of partitions
  for (n in seq_along(partition)){
    
    # Set the required values for this section
    part <- partition[[n]]
    idPart <- which(chromosomes %in% part)
    chrEnd <- lengthSeq[chromosomes %in% part]
    
    # Get & Summarise the reads
    # Return the reads from the bam file as a list, with each element containing reads from a single chr
    sbp <- ScanBamParam(what = scanWhat,which = GRanges(seqnames = part,ranges = IRanges(start = 1,end = chrEnd)))
    readInfo <- scanBam(file, param = sbp)
    if (mapqFilter > 0){
      #get only reads that pass the mapping quality filter
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
    
    if (sum(nReadsInPart) > 0){
      # Calculate the windows that overlap mustKeepRanges
      mustKeepWin <- list()
      if (!missing(mustKeepRanges)){
        if (length(intersect(allChromosomesMustKeep,part))>0){
          mustKeepWin <- getWinFromGranges(mustKeepRanges[seqnames(mustKeepRanges) %in% part],chromosomeInfo[idPart,],winWidth,winStep)
        }
      }
      
      # Calculate the first/last bases/reads in the chromosome partition
      chromosomeInfo[idPart,] <- chromosomeInfoInPartition(chromosomeInfo[idPart,], winStep)
      
      # Concatenate several lists of the chromosome partition into one list
      readInfo <- concatenateAlignments(readInfo,chromosomeInfo[idPart,])
      
      # Get the index of R1 and R2 reads and process each subset separately
      if (paired){
        firstReadIndex <- ((floor(readInfo$flag/64) %% 2) == 1)
        secondReadIndex <- !firstReadIndex
        subset <- list("R1"=firstReadIndex,"R2"=secondReadIndex)
      } else{
        subset <- list(NULL)
      }
      for (s in seq_along(subset)){
        winPositiveAlignments <- getWinOfAlignments(readInfo,"+",winWidth,winStep,readProp = readProp, useCoverage= (getWin || useCoverage),subset[[s]])
        winNegativeAlignments <- getWinOfAlignments(readInfo,"-",winWidth,winStep,readProp = readProp, useCoverage= (getWin || useCoverage),subset[[s]])
        
        probaWin <- keptProbaWin(winPositiveAlignments,winNegativeAlignments,winWidth,winStep,logitThreshold,pvalueThreshold,errorRate,mustKeepWin,min,max,getWin,useCoverage=useCoverage) # the probability of each positive/negative window; this probability will be assigned to every positive/negative read in that window
        if (getWin){
          win <- getWinInChromosome(probaWin$Win,part,chromosomeInfo[idPart,],winWidth,winStep)
          if (s==1){
            allWin[[n]] <- win
          } else{
            allWin[[n]]$Type <- Rle("R1")
            win$Type <- Rle("R2")
            allWin[[n]] <- rbind(allWin[[n]],win)
          }
        }
        
        # the positive read fragments to be kept
        keptPositiveAlignment <- keptReadFragment(winPositiveAlignments$Win,probaWin$Positive,errorRate) 
        # the negative read fragments to be kept
        keptNegativeAlignment <- keptReadFragment(winNegativeAlignments$Win,probaWin$Negative,errorRate) 
       
        rm(probaWin)
        
        # Infer the index of kept positive/negative alignments 
        if (s==1){
          keptAlignments <- c(unique(mcols(winPositiveAlignments$Win)$alignment[keptPositiveAlignment]),unique(mcols(winNegativeAlignments$Win)$alignment[keptNegativeAlignment]))  # the vector of all alignments to be kept  
        } else{
          kept <- c(unique(mcols(winPositiveAlignments$Win)$alignment[keptPositiveAlignment]),unique(mcols(winNegativeAlignments$Win)$alignment[keptNegativeAlignment]))  # the vector of all alignments to be kept
          keptAlignments <- sort(c(keptAlignments,kept))
        }
        rm(winPositiveAlignments,winNegativeAlignments)
      }
      rm(readInfo)
      
      #Write the kept alignments into the output bamfile
      #Run through each chromsome
      for (i in seq_along(part)){
        id <- which(chromosomes == part[i])
        if (chromosomeInfo$NbOriginalReads[id]>0){
          #Get the kept alignments of the considering chromosome 
          rangeInChr <- which((keptAlignments>=chromosomeInfo$FirstReadInPartition[id] & keptAlignments<= chromosomeInfo$LastReadInPartition[id]))
          chromosomeInfo$NbKeptReads[id] <- length(rangeInChr)
          if (chromosomeInfo$NbKeptReads[id]>0){
            index <- keptAlignments[rangeInChr]- chromosomeInfo$FirstReadInPartition[id]+1
            if (mapqFilter>0) index <- mqfilter[[i]][index]
            #write the kept alignments into the output bam file  
            bamSave(writer,bamRange(reader,c(which(allRefChromosomes == part[i])-1,0,lengthSeq[id]))[index,])
          }
          rm(rangeInChr)
        }
      }
      cat(paste0("Sequence ",part,", length: ",chromosomeInfo$Length[idPart],
                 ", number of reads: ",chromosomeInfo$NbOriginalReads[idPart],
                 ", number of kept reads: ",chromosomeInfo$NbKeptReads[idPart],"\n"),file=statfile,append=append)
      if (!append) append <- TRUE
      rm(keptAlignments)
    }
  }
  bamClose(writer)
  bamClose(reader)
  
  cat("Summary:\n",file = statfile, append = append)
  cat(paste0("Number of original reads: ",sum(chromosomeInfo$NbOriginalReads),", number of kept reads: ",sum(chromosomeInfo$NbKeptReads),", removal proportion: ",(sum(chromosomeInfo$NbOriginalReads)-sum(chromosomeInfo$NbKeptReads))/sum(chromosomeInfo$NbOriginalReads)),"\n",file = statfile,append=TRUE)
  endTime <- proc.time()
  cat(paste0("Total elapsed time ",(endTime-startTime)[[3]]/60," minutes\n"),file = statfile,append=TRUE)
  if (getWin){
    return(do.call(rbind,allWin))
  }
}



