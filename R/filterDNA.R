#' @title Filter Single End Bam File
#'
#' @description Filter putative double strand DNA from a strand specific single-end RNA-seq using a window sliding across the genome.

#'
#' @param bamfilein the input single end bam file to be filterd. Your bamfile should be sorted and have an index file located at the same path as well.
#' @param bamfileout the output filtered bam file
#' @param statfile the file to write the summary of the results
#' @param chromosomes the list of chromosomes to be filtered
#' @param yieldSize by default is 1e8, i.e. the bam file is read by block of chromosomes such that the total length of each block is at least 1e8
#' @param mustKeepRanges a GRanges object defines the ranges such that every read maps to those ranges must be always kept regardless the strand proportion of the windows containing them.
#' @param getWin if TRUE, the function will return a data frame containing the information of all windows. It's FALSE by default.
#' @param win the length of the sliding window, 1000 by default.
#' @param step the step length to sliding the window, 100 by default.
#' @param threshold the threshold upper which we keep the windows. 0.7 by default
#' @param pvalueThreshold the threshold for the p-value. 0.05 by default
#' @param min if a window has least than \code{min} reads, then it will be rejected regardless the strand proportion. 0 by default
#' @param max if a window has more than \code{max} reads, then it will be kept regardless the strand proportion. If 0 then it doesn't have effect on selecting window. 0 by default.
#' @param errorRate the probability that an RNA read takes the false strand. 0.01 by default
#' @param limit a read is considered to be included in a window if and only if at least \code{limit} percent of it is in the window. 0.75 by default
#' @param coverage if TRUE, then the strand information in each window corresponds to the sum of coverage coming from positive/negative reads; and not the number of positive/negative reads as default.
#'
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
#' @seealso filterDNAPairs, getWinFromBamFile, getWinFromPairedBamFile, plotHist, plotWin
#'
#' @examples
#' bamfilein <- system.file("data","s1.chr1.bam",package = "rnaCleanR")
#' filterDNA(bamfilein,bamfileout="out.bam",statfile = "out.stat")
#' @export
#' @importFrom rbamtools bamReader getHeader bamWriter bamClose bamSave bamRange
#' @importFrom Rsamtools BamFile scanBam ScanBamParam
#' @importFrom GenomicRanges GRanges ranges
#' @importFrom GenomeInfoDb seqinfo seqnames seqlengths
#' @importFrom IRanges IRanges
#' @importFrom magrittr %>%
#' @import S4Vectors
#'
filterDNA <- function(bamfilein,bamfileout,statfile,chromosomes,yieldSize = 1e8, mustKeepRanges,getWin=FALSE,win=1000,step=100,threshold=0.7,pvalueThreshold=0.05,min=0,max=0,errorRate=0.01,limit=0.75,coverage=FALSE){
  startTime <- proc.time()
  bf <- BamFile(bamfilein)
  seqinfo <- seqinfo(bf)
  allChromosomes <- seqnames(seqinfo)
  lengthSeq <- seqlengths(seqinfo)
  if (missing(chromosomes)){
    chromosomes <- allChromosomes
  }
  partition <- partitionChromosomes(chromosomes,lengthSeq[allChromosomes %in% chromosomes],yieldSize = yieldSize)

  reader <- bamReader(bamfilein,idx=TRUE) #open a reader of the input bamfile to extract read afterward
  header <- getHeader(reader) #get the header of the input bam file
  writer <- bamWriter(header,bamfileout) #prepare to write the output bamfile with the same header
  logitThreshold <- binomial()$linkfun(threshold)

  if (missing(statfile)){
    message("Summary will be written to file out.stat")
    statfile <- "out.stat"
  }
  
  if (!missing(mustKeepRanges)){
    allChromosomesMustKeep <- levels(seqnames(mustKeepRanges))
  }
  statInfo <- data.frame("Sequence"="chr","Length"=rep(0,length(chromosomes)),
                         "NbOriginalReads" = rep(0,length(chromosomes)), 
                         "NbKeptReads" = rep(0,length(chromosomes)),
                         "FirstBaseInPartition" = rep(NA,length(chromosomes)),
                         "LastBaseInPartition" = rep(NA,length(chromosomes)),
                         "FirstReadInPartition" = rep(NA,length(chromosomes)),
                         "LastReadInPartition" = rep(NA,length(chromosomes)),
                         stringsAsFactors = FALSE)

  if (coverage==TRUE){
    allWin <- data.frame("Chr"=c(), "Start" = c(), "NbPositive"= c(), "NbNegative"= c(),"MaxCoverage" = c()) 
  }
  else{
    allWin <- data.frame("Chr"=c(), "Start" = c(), "NbPositive"= c(), "NbNegative"= c())  
  }
  append <- FALSE
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
      mustKeepWin <- list()
      if (!missing(mustKeepRanges)){
        if (length(intersect(allChromosomesMustKeep,part))>0){
          mustKeepWin <- getWinFromGranges(mustKeepRanges[seqnames(mustKeepRanges) %in% part],part,statInfo[idPart,],win,step)
        }
      }
      
      winPositiveAlignments <- getWinOfAlignments(bam,"+",win,step,limit,coverage=coverage)
      winNegativeAlignments <- getWinOfAlignments(bam,"-",win,step,limit,coverage=coverage)
      rm(bam)

      probaWin <- keptProbaWin(winPositiveAlignments,winNegativeAlignments,win,step,logitThreshold,pvalueThreshold,errorRate,mustKeepWin,min,max,getWin,coverage=coverage) # the probability of each positive/negative window; this probability will be assigned to every positive/negative read in that window
      keptPositiveAlignment <- keptAlignment(winPositiveAlignments$Win,probaWin$Positive,errorRate) # the positive alignments to be kept
      keptNegativeAlignment <- keptAlignment(winNegativeAlignments$Win,probaWin$Negative,errorRate) # the negative alignments to be kept

      if (getWin){
        allWin <- rbind(allWin,getWinInChromosome(probaWin$Win,part,statInfo[idPart,],win,step))
      }
      rm(probaWin)
      
      keptAlignment <- c(unique(mcols(winPositiveAlignments$Win)$alignment[keptPositiveAlignment]),unique(mcols(winNegativeAlignments$Win)$alignment[keptNegativeAlignment])) %>% sort() # the vector of all alignments to be kept
      rm(winPositiveAlignments,winNegativeAlignments)

      for (i in seq_along(part)){
        id <- which(chromosomes == part[i])
        if (statInfo$NbOriginalReads[id]>0){
          rangeInChr <- which((keptAlignment>=statInfo$FirstReadInPartition[id] & keptAlignment<= statInfo$LastReadInPartition[id]))
          statInfo$NbKeptReads[id] <- length(rangeInChr)
        }
        if (statInfo$NbKeptReads[id]>0){
          chromosomeIndex <- which(allChromosomes == part[i])
          bamSave(writer,bamRange(reader,c(chromosomeIndex-1,0,lengthSeq[chromosomeIndex]))[keptAlignment[rangeInChr]- statInfo$FirstReadInPartition[id]+1,],refid=chromosomeIndex-1)#write the kept alignments into the output bam file
          rm(rangeInChr)
        }
      }
      cat(paste0("Sequence ",part,", length: ",statInfo$Length[idPart],
                 ", number of reads: ",statInfo$NbOriginalReads[idPart],
                 ", number of kept reads: ",statInfo$NbKeptReads[idPart],"\n"),file=statfile,append=append)
      if (!append) append <- TRUE
      rm(keptAlignment)
    }
  }
  bamClose(writer)
  bamClose(reader)

  cat("Summary:\n",file = statfile, append = append)
  cat(paste0("Number of original reads: ",sum(statInfo$NbOriginalReads),", number of kept reads: ",sum(statInfo$NbKeptReads),", removal proportion: ",(sum(statInfo$NbOriginalReads)-sum(statInfo$NbKeptReads))/sum(statInfo$NbOriginalReads)),"\n",file = statfile,append=TRUE)
  endTime <- proc.time()
  cat(paste0("Total elapsed time ",(endTime-startTime)[[3]]/60," minutes\n"),file = statfile,append=TRUE)
  if (getWin){
    return(allWin %>% dplyr::mutate("Start" = (Start-1)*step+1))
  }
}



