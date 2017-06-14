#' @title Filter Single End Bam File
#'
#' @description Filter putative double strand DNA from a strand specific RNA-seq using a window sliding across the genome.

#'
#' @param bamfilein the input bam file to be filterd. Your bamfile should be sorted and have an index file located at the same path as well.
#' @param bamfileout the output filtered bam file
#' @param statfile the file to write the summary of the results
#' @param chromosomes the list of chromosomes to be filtered
#' @param yieldSize by default is 1e8, i.e. the bam file is read by block of chromosomes such that the total length of each block is at least 1e8
#' @param mustKeepRanges a GRanges object defines the ranges such that every read maps to those ranges must be always kept regardless the strand proportion of the windows containing them.
#' @param getWin if TRUE, the function will return a data frame containing the information of all windows. It's FALSE by default.
#' @param win the length of the sliding window, 1000 by default.
#' @param step the step length to sliding the window, 100 by default.
#' @param threshold the threshold upper which we keep the reads. 0.7 by default
#' @param pvalueThreshold the threshold for the p-value. 0.05 by default
#' @param min if a window has least than min reads, then it will be rejected regardless the strand proportion. 0 by default
#' @param max if a window has more than max reads, then it will be kept regardless the strand proportion. If 0 then it doesn't have effect on selecting window. 0 by default.
#' @param errorRate the probability that an RNA read takes the false strand. 0.01 by default
#' @param limit a read is considered to be included in a window if and only if at least limit percent of it is in the window. 0.75 by default
#' @param coverage if TRUE, then the strand information in each window corresponds to the sum of coverage coming from positive/negative reads; and not the number of positive/negative reads as default.
#'
#'
#' @details filterDNA reads a bam file containing strand specific RNA reads, and filter putative double strand DNA.
#' Using a window sliding across the genome, we calculate the positive/negative proportion of reads in that window.
#' For each window, we use logistic regression to estimate the proportion of reads in the window derived from
#' stranded RNA (positive or negative).
#'
#' \eqn{\pi}: proportion of reads in the window derived from stranded RNA (positive or negative)
#'
#' Null hypothesis: \eqn{\pi < {\pi}_{0}} where \eqn{{\pi}_{0}} is the given threshold.
#'
#' Only windows with p-value <= 0.05 are kept. Considering a positive window that is kept, let P be its number of positive reads, and let M
#' be its number of negative reads. Since these M negative reads should come from double-strand DNA, then there should be also M postive reads among the
#' P positive reads come from double-strand DNA. In other words, there are only (P-M) positive reads come from RNA. Hence, each positive read in this window is kept
#' with the probability equalling (P-M)/P. Each negative read is kept with the probability equalling the given errorRate which is the rate that
#' an RNA read of your sample has wrong strand.
#'
#' Since each alignment can be belonged to several windows, then the probability of keeping an alignment is the maximum probability defined by
#' all windows that contain it.
#'
#' @seealso filterDNAPairs, getWin, getWinPairs, plotHist, plotWin
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

  nbOriginalReads <- 0 #number of original reads
  nbKeptReads <- 0 #number of kept reads

  mustKeepWin <- list()
  if (!missing(mustKeepRanges)){
    mustKeepWin <- getWinFromGranges(mustKeepRanges,chromosomes,win,step) # the windows must be kept, calculated from the input must be kept ranges
  }

  allWin <- data.frame("Chr"=c(), "Start" = c(), "NbPositiveReads"= c(), "NbNegativeReads"= c())
  append <- FALSE
  for (part in partition){
    lengthSeqInChr <- lengthSeq[allChromosomes %in% part]
    lengthSeqInPart <- c(0,cumsum(as.numeric(lengthSeqInChr)))
    lengthSeqInPart <- step*ceiling(lengthSeqInPart/step)

    bam <- scanBam(bamfilein,param=ScanBamParam(what=c("pos","cigar","strand"),
                                                which=GRanges(seqnames = part,ranges = IRanges(start=1,end=lengthSeqInChr))))

    nbOriginalReadsInChr <- sapply(1:length(bam),function(i){length(bam[[i]]$strand)})
    if (sum(nbOriginalReadsInChr)>0){
      nil <- which(nbOriginalReadsInChr!=0)
      nbOriginalReadsInChr <- nbOriginalReadsInChr[nil]
      part <- part[nil]
      nbOriginalReadsInPart <- c(0,cumsum(nbOriginalReadsInChr))
      lengthSeqInPart <- lengthSeqInPart[c(1,nil+1)]
      bam <- concatenateAlignments(bam[nil],nbOriginalReadsInPart,lengthSeqInPart,sum(nbOriginalReadsInChr))
      winPositiveAlignments <- getWinOfAlignments(bam,"+",win,step,limit,coverage=coverage)
      winNegativeAlignments <- getWinOfAlignments(bam,"-",win,step,limit,coverage=coverage)
      rm(bam)

      probaWin <- keptProbaWin(winPositiveAlignments,winNegativeAlignments,logitThreshold,errorRate,mustKeepWin,min,max,getWin,coverage=coverage) # the probability of each positive/negative window; this probability will be assigned to every positive/negative read in that window
      keptPositiveAlignment <- keptAlignment(winPositiveAlignments$Win,probaWin$Positive,errorRate) # the positive alignments to be kept
      keptNegativeAlignment <- keptAlignment(winNegativeAlignments$Win,probaWin$Negative,errorRate) # the negative alignments to be kept

      if (getWin){
        Chromosome <- rep("",nrow(probaWin$Win))
        for (i in seq_along(part)){
          mi <- ceiling((lengthSeqInPart[i]+1)/step)
          ma <- ceiling((lengthSeqInPart[i+1]-win+1)/step)
          j <- which(probaWin$Win$Start >=mi & probaWin$Win$Start <=ma)
          Chromosome[j] <- part[i]
          probaWin$Win$Start[j] <- probaWin$Win$Start[j] - mi +1
        }
        probaWin$Win[["Chr"]] <- Chromosome
        allWin <- rbind(allWin,probaWin$Win)
      }
      rm(probaWin)
      keptAlignment <- c(unique(mcols(winPositiveAlignments$Win)$alignment[keptPositiveAlignment]),unique(mcols(winNegativeAlignments$Win)$alignment[keptNegativeAlignment])) %>% sort() # the vector of all alignments to be kept
      rm(winPositiveAlignments,winNegativeAlignments)

      nbKeptReadsInChr <- rep(0,length(nbOriginalReadsInChr))
      for (i in 1:length(nbOriginalReadsInChr)){
        chromosomeIndex <- which(allChromosomes == part[i])
        range <- (keptAlignment>nbOriginalReadsInPart[i] & keptAlignment<=nbOriginalReadsInPart[i+1])
        nbKeptReadsInChr[i] <- sum(range)
        if (nbKeptReadsInChr[i]>0){
          bamSave(writer,bamRange(reader,c(chromosomeIndex-1,0,lengthSeq[chromosomeIndex]))[keptAlignment[range]-nbOriginalReadsInPart[i],],refid=chromosomeIndex-1)#write the kept alignments into the output bam file
        }
        rm(range)
        cat(paste0("Sequence ",part[i],", length: ",lengthSeq[chromosomeIndex],", number of reads: ",nbOriginalReadsInChr[i],", number of kept reads: ",nbKeptReadsInChr[i],"\n"),file=statfile,append=append)
        if (!append) append <- TRUE
      }
      rm(keptAlignment)
      nbOriginalReads <- nbOriginalReads + sum(nbOriginalReadsInChr)
      nbKeptReads <- nbKeptReads + sum(nbKeptReadsInChr)
    }
  }
  bamClose(writer)
  bamClose(reader)

  cat("Summary:\n",file = statfile, append = append)
  cat(paste0("Number of original reads: ",nbOriginalReads,", number of kept reads: ",nbKeptReads,", removal proportion: ",(nbOriginalReads-nbKeptReads)/nbOriginalReads),"\n",file = statfile,append=TRUE)
  endTime <- proc.time()
  cat(paste0("Total elapsed time ",(endTime-startTime)[[3]]/60," minutes\n"),file = statfile,append=TRUE)
  if (getWin){
    return(allWin %>% dplyr::mutate("Start" = (Start-1)*step+1))
  }
}



