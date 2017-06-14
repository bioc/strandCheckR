#' @title Filter Paired End Bam File
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
#' @param minReads if a window has least than minReads reads, then it will be rejected regardless the strand proportion. 0 by default
#' @param maxReads if a window has more than maxReads reads, then it will be kept regardless the strand proportion. If 0 then it doesn't have effect on selecting window. 0 by default.
#' @param errorRate the probability that an RNA read takes the false strand. 0.01 by default
#' @param limit a read is considered to be included in a window if and only if at least limit percent of it is in the window. 0.75 by default
#' @param pair if "free" than the any pair of reads do not need to be both kept in the filtered file, i.e. first reads and second reads are filtered independently. Otherwise, if pair = "intersection" then a pair of reads is kept if and only if both first and second reads are kept; if pair = "union" then a pair of reads is kept if either the first or the second read is kept.
#' @param coverage if TRUE, then the strand information in each window corresponds to the sum of coverage coming from positive/negative reads; and not the number of positive/negative reads as default.
#'
#' @details filterDNA reads a bam file containing strand specific RNA reads, and filter putative double strand DNA.
#' Using a window sliding across the genome, we calculate the positive/negative proportion of reads in that window.
#' For each window, we use logistic regression to estimate the proportion of reads in the window derived from
#' stranded RNA (positive or negative).
#' @seealso filterDNA, getWin, getWinPairs, plotHist, plotWin
#'
#' @examples
#' bamfilein <- system.file("data","120.10.bam",package = "rnaCleanR")
#' filterDNAPairs(bamfilein,bamfileout="out.bam",statfile = "out.stat")
#' @export
#' @importFrom rbamtools bamReader getHeader bamWriter bamClose bamSave bamRange
#' @importFrom Rsamtools BamFile scanBam ScanBamParam
#' @importFrom GenomicRanges GRanges ranges
#' @importFrom GenomeInfoDb seqinfo seqnames seqlengths
#' @importFrom IRanges IRanges
#' @importFrom magrittr %>%
#' @import S4Vectors
#'
filterDNAPairs <- function(bamfilein,bamfileout,statfile,chromosomes,yieldSize = 1e8,mustKeepRanges,getWin=FALSE,win=1000,step=100,threshold=0.7,pvalueThreshold=0.05,minReads=0,maxReads=0,errorRate=0.01,limit=0.75,pair="free",coverage=FALSE){
  stopifnot(pair=="free" | pair=="intersection" | pair=="union")
  startTime <- proc.time()
  bf <- BamFile(bamfilein)
  seqinfo <- seqinfo(bf)
  allChromosomes <- seqnames(seqinfo)
  lengthSeq <- seqlengths(seqinfo)
  if (missing(chromosomes)){
    chromosomes <- allChromosomes
  }
  partition <- partitionChromosomes(chromosomes,lengthSeq[allChromosomes %in% chromosomes],yieldSize)
  reader <- bamReader(bamfilein,idx=TRUE) #open a reader of the input bamfile to extract read afterward
  header <- getHeader(reader) #get the header of the input bam file
  writer <- bamWriter(header,bamfileout) #prepare to write the output bamfile with the same header
  logitThreshold <- binomial()$linkfun(threshold)

  if (missing(statfile)){
    message("Summary will be written to file out.stat")
    statfile <- "out.stat"
  }

  nbOriginalFirstReads <- 0 #number of original reads
  nbOriginalSecondReads <- 0
  nbKeptFirstReads <- 0 #number of kept reads
  nbKeptSecondReads <- 0
  mustKeepWin <- list()
  if (!missing(mustKeepRanges)){
    mustKeepWin <- getWinFromGranges(mustKeepRanges,chromosomes,win,step) # the windows must be kept, calculated from the input must be kept ranges
  }
  allWin <- data.frame("Type"=c(),"Chr"=c(), "Start" = c(), "NbPositiveReads"= c(), "NbNegativeReads"= c())
  append <- FALSE
  for (part in partition){
    lengthSeqInChr <- lengthSeq[allChromosomes %in% part]
    lengthSeqInPart <- c(0,cumsum(as.numeric(lengthSeqInChr)))
    lengthSeqInPart <- step*ceiling(lengthSeqInPart/step)

    if (pair=="free"){
      bam <- scanBam(bamfilein,param=ScanBamParam(what=c("pos","cigar","strand","flag"),
                                                  which=GRanges(seqnames = part,ranges = IRanges(start=1,end=lengthSeq[allChromosomes %in% part]))))
    }
    else{
      bam <- scanBam(bamfilein,param=ScanBamParam(what=c("pos","cigar","strand","flag","qname"),
                                                  which=GRanges(seqnames = part,ranges = IRanges(start=1,end=lengthSeq[allChromosomes %in% part]))))
    }

    nbOriginalReadsInChr <- sapply(1:length(bam),function(i){length(bam[[i]]$strand)})
    if (sum(nbOriginalReadsInChr)>0){
      nil <- which(nbOriginalReadsInChr!=0)
      nbOriginalReadsInChr <- nbOriginalReadsInChr[nil]
      part <- part[nil]
      nbOriginalReadsInPart <- c(0,cumsum(nbOriginalReadsInChr))
      lengthSeqInPart <- lengthSeqInPart[c(1,nil+1)]
      if (pair=="free"){
        bam <- concatenateAlignments(bam[nil],nbOriginalReadsInPart,lengthSeqInPart,sum(nbOriginalReadsInChr),flag=TRUE)
      }
      else{
        bam <- concatenateAlignments(bam[nil],nbOriginalReadsInPart,lengthSeqInPart,sum(nbOriginalReadsInChr),flag=TRUE,qname=TRUE)
      }

      firstReadIndex <- ((floor(bam$flag/64) %% 2) == 1)
      secondReadIndex <- !firstReadIndex
      winFirstPositiveAlignments <- getWinOfAlignments(bam,"+",win,step,limit,firstReadIndex,coverage=coverage)
      winSecondPositiveAlignments <- getWinOfAlignments(bam,"+",win,step,limit,secondReadIndex,coverage=coverage)
      winFirstNegativeAlignments <- getWinOfAlignments(bam,"-",win,step,limit,firstReadIndex,coverage=coverage)
      winSecondNegativeAlignments <- getWinOfAlignments(bam,"-",win,step,limit,secondReadIndex,coverage=coverage)

      probaWinFirst <- keptProbaWin(winFirstPositiveAlignments,winFirstNegativeAlignments,logitThreshold,errorRate,mustKeepWin,minReads,maxReads,getWin,coverage=coverage) # the probability of each positive/negative window; this probability will be assigned to every positive/negative read in that window
      probaWinSecond <- keptProbaWin(winSecondPositiveAlignments,winSecondNegativeAlignments,logitThreshold,errorRate,mustKeepWin,minReads,maxReads,getWin,coverage=coverage) # the probability of each positive/negative window; this probability will be assigned to every positive/negative read in that window
      if (getWin){
        ChromosomeFirst <- rep("",nrow(probaWinFirst$Win))
        ChromosomeSecond <- rep("",nrow(probaWinSecond$Win))
        for (i in seq_along(part)){
          mi <- ceiling((lengthSeqInPart[i]+1)/step)
          ma <- ceiling((lengthSeqInPart[i+1]-win+1)/step)

          j <- which(probaWinFirst$Win$Start >=mi & probaWinFirst$Win$Start <=ma)
          ChromosomeFirst[j] <- part[i]
          probaWinFirst$Win$Start[j] <- probaWinFirst$Win$Start[j] - mi +1

          j <- which(probaWinSecond$Win$Start >=mi & probaWinSecond$Win$Start <=ma)
          ChromosomeSecond[j] <- part[i]
          probaWinSecond$Win$Start[j] <- probaWinSecond$Win$Start[j] - mi +1
        }
        probaWinFirst$Win[["Chr"]] <- ChromosomeFirst
        probaWinSecond$Win[["Chr"]] <- ChromosomeSecond
        allWin <- rbind(allWin,data.frame("Type"="First",probaWinFirst$Win),data.frame("Type"="Second",probaWinSecond$Win))
      }
      keptFirstPositiveAlignment <- keptAlignment(winFirstPositiveAlignments$Win,probaWinFirst$Positive,errorRate) # the positive alignments to be kept
      keptFirstNegativeAlignment <- keptAlignment(winFirstNegativeAlignments$Win,probaWinFirst$Negative,errorRate) # the negative alignments to be kept
      keptSecondPositiveAlignment <- keptAlignment(winSecondPositiveAlignments$Win,probaWinSecond$Positive,errorRate) # the positive alignments to be kept
      keptSecondNegativeAlignment <- keptAlignment(winSecondNegativeAlignments$Win,probaWinSecond$Negative,errorRate) # the negative alignments to be kept
      rm(probaWinFirst,probaWinSecond)
      keptFirstAlignments <- c(unique(mcols(winFirstPositiveAlignments$Win)$alignment[keptFirstPositiveAlignment]),unique(mcols(winFirstNegativeAlignments$Win)$alignment[keptFirstNegativeAlignment]))  # the vector of all alignments to be kept
      keptSecondAlignments <- c(unique(mcols(winSecondPositiveAlignments$Win)$alignment[keptSecondPositiveAlignment]),unique(mcols(winSecondNegativeAlignments$Win)$alignment[keptSecondNegativeAlignment])) # the vector of all alignments to be kept
      rm(winFirstPositiveAlignments,winFirstNegativeAlignments,winSecondPositiveAlignments,winSecondNegativeAlignments)

      if (pair=="free"){
        keptAlignments <- sort(c(keptFirstAlignments,keptSecondAlignments))
      }
      else{
        if (pair=="intersection"){
          keptReadNames <- intersect(bam$qname[keptFirstAlignments],bam$qname[keptSecondAlignments])
        }
        else if (pair=="union"){
          commonNames <- intersect(bam$qname[firstReadIndex],bam$qname[secondReadIndex])
          keptReadNames <-intersect(commonNames,union(bam$qname[keptFirstAlignments],bam$qname[keptSecondAlignments]))
        }
        keptAlignments <- which(bam$qname %in% keptReadNames)
      }
      rm(bam)
      nbOriginalFirstReadsInChr <- rep(0,length(nbOriginalReadsInChr))
      nbOriginalSecondReadsInChr <- rep(0,length(nbOriginalReadsInChr))
      nbKeptFirstReadsInChr <- rep(0,length(nbOriginalReadsInChr))
      nbKeptSecondReadsInChr <- rep(0,length(nbOriginalReadsInChr))
      for (i in 1:length(nbOriginalReadsInChr)){
        range <- keptAlignments[(keptAlignments>nbOriginalReadsInPart[i] & keptAlignments<=nbOriginalReadsInPart[i+1])]
        nbKeptFirstReadsInChr[i] <- sum(firstReadIndex[range])
        nbKeptSecondReadsInChr[i] <- sum(secondReadIndex[range])
        chromosomeIndex <- which(allChromosomes == part[i])
        if (nbKeptFirstReadsInChr[i]>0 | nbKeptSecondReadsInChr[i]>0){
          nbOriginalFirstReadsInChr[i] <- sum(firstReadIndex[(1+nbOriginalReadsInPart[i]):nbOriginalReadsInPart[i+1]])
          nbOriginalSecondReadsInChr[i] <- sum(secondReadIndex[(1+nbOriginalReadsInPart[i]):nbOriginalReadsInPart[i+1]])
          bamSave(writer,bamRange(reader,c(chromosomeIndex-1,0,lengthSeq[chromosomeIndex]))[range-nbOriginalReadsInPart[i],],refid=chromosomeIndex-1)#write the kept alignments into the output bam file
        }
        rm(range)
        cat(paste0("Sequence ",part[i],", length: ",lengthSeq[chromosomeIndex],", number of first reads: ",nbOriginalFirstReadsInChr[i],", number of second reads: ",nbOriginalSecondReadsInChr[i],", number of kept first reads: ",nbKeptFirstReadsInChr[i],", number of kept second reads: ",nbKeptSecondReadsInChr[i],"\n"),file=statfile,append=append)
        if (!append) append <- TRUE
      }
      rm(keptAlignments)
      nbOriginalFirstReads <- nbOriginalFirstReads + sum(nbOriginalFirstReadsInChr)
      nbOriginalSecondReads <- nbOriginalSecondReads + sum(nbOriginalSecondReadsInChr)
      nbKeptFirstReads <- nbKeptFirstReads + sum(nbKeptFirstReadsInChr)
      nbKeptSecondReads <- nbKeptSecondReads + sum(nbKeptSecondReadsInChr)
    }
  }
  bamClose(writer)
  bamClose(reader)
  cat("Summary:\n",file = statfile, append = append)
  cat(paste0("Number of original first reads: ",nbOriginalFirstReads,", number of original second reads: ",nbOriginalSecondReads,", number of kept first reads: ",nbKeptFirstReads,", number of kept second reads: ",nbKeptSecondReads,", removal proportion of first reads: ",(nbOriginalFirstReads-nbKeptFirstReads)/nbOriginalFirstReads,", removal proportion of second reads: ",(nbOriginalSecondReads-nbKeptSecondReads)/nbOriginalSecondReads),"\n",file = statfile,append=TRUE)
  endTime <- proc.time()
  cat(paste0("Total elapsed time ",(endTime-startTime)[[3]]/60," minutes\n"),file = statfile,append=TRUE)
  if (getWin){
    return(allWin %>% dplyr::mutate("Start" = (Start-1)*step+1))
  }
}



