#' @title Filter Paired End Bam File
#'
#' @description Filter putative double strand DNA from a strand specific paried-end RNA-seq using a window sliding across the genome.

#'
#' @param bamfilein the input paired-end bam file to be filterd. Your bamfile should be sorted and have an index file located at the same path as well.
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
#' @param min if a window has least than \code{min} reads, then it will be rejected regardless the strand proportion. 0 by default
#' @param max if a window has more than \code{max} reads, then it will be kept regardless the strand proportion. If 0 then it doesn't have effect on selecting window. 0 by default.
#' @param errorRate the probability that an RNA read takes the false strand. 0.01 by default
#' @param limit a read is considered to be included in a window if and only if at least \code{limit} percent of it is in the window. 0.75 by default
#' @param pair if "free" than the any pair of reads do not need to be both kept in the filtered file, i.e. first reads and second reads are filtered independently. Otherwise, if pair = "intersection" then a pair of reads is kept if and only if both first and second reads are kept; if pair = "union" then a pair of reads is kept if either the first or the second read is kept.
#' @param coverage if TRUE, then the strand information in each window corresponds to the sum of coverage coming from positive/negative reads; and not the number of positive/negative reads as default.
#'
#' @details filterDNAPairs reads a paired-end bam file containing strand specific paired-end RNA reads, and filter putative double strand DNA.
#' Using a window sliding across the genome, we calculate the positive/negative proportion of reads in that window.
#' For each window, we use logistic regression to estimate the proportion of reads in the window derived from
#' stranded RNA (positive or negative).
#' @seealso filterDNA, getWinFromBamFile, getWinFromPairedBamFile, plotHist, plotWin
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
filterDNAPairs <- function(bamfilein,bamfileout,statfile,chromosomes,yieldSize = 1e8,mustKeepRanges,getWin=FALSE,win=1000,step=100,threshold=0.7,pvalueThreshold=0.05,min=0,max=0,errorRate=0.01,limit=0.75,pair="free",coverage=FALSE){
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
  if (!missing(mustKeepRanges)){
    allChromosomesMustKeep <- levels(seqnames(mustKeepRanges))
  }
  statInfo <- data.frame("Sequence"="chr","Length"=rep(0,length(chromosomes)),
                         "NbOriginalReads" = rep(0,length(chromosomes)), 
                         "NbOriginalFirstReads" = rep(0,length(chromosomes)), 
                         "NbOriginalSecondReads" = rep(0,length(chromosomes)), 
                         "NbKeptFirstReads" = rep(0,length(chromosomes)),
                         "NbKeptSecondReads" = rep(0,length(chromosomes)),
                         "FirstBaseInPartition" = rep(NA,length(chromosomes)),
                         "LastBaseInPartition" = rep(NA,length(chromosomes)),
                         "FirstReadInPartition" = rep(NA,length(chromosomes)),
                         "LastReadInPartition" = rep(NA,length(chromosomes)),
                         stringsAsFactors = FALSE)
 if (coverage){
   allWin <- data.frame("Type"=c(),"Chr"=c(), "Start" = c(), "NbPositive"= c(), "NbNegative"= c(), "MaxCoverage" = c())
 } else{
   allWin <- data.frame("Type"=c(),"Chr"=c(), "Start" = c(), "NbPositive"= c(), "NbNegative"= c())
 }
  append <- FALSE
  for (part in partition){
    idPart <- which(chromosomes %in% part)
    statInfo$Sequence[idPart] <- part
    statInfo$Length[idPart] <- lengthSeq[allChromosomes %in% part]

    if (pair=="free"){
      bam <- scanBam(bamfilein,param=ScanBamParam(what=c("pos","cigar","strand","flag"),
                                                  which=GRanges(seqnames = part,ranges = IRanges(start=1,end=statInfo$Length[idPart]))))
    }
    else{
      bam <- scanBam(bamfilein,param=ScanBamParam(what=c("pos","cigar","strand","flag","qname"),
                                                  which=GRanges(seqnames = part,ranges = IRanges(start=1,end=statInfo$Length[idPart]))))
    }

    statInfo$NbOriginalReads[idPart] <- sapply(seq_along(bam),function(i){length(bam[[i]]$strand)})
    if (sum(statInfo$NbOriginalReads[idPart])>0){
      statInfo[idPart,] <- statInfoInPartition(statInfo[idPart,],step)
      if (pair=="free"){
        bam <- concatenateAlignments(bam,statInfo[idPart,],flag=TRUE)
      }
      else{
        bam <- concatenateAlignments(bam,statInfo[idPart,],flag=TRUE,qname=TRUE)
      }
      mustKeepWin <- list()
      if (!missing(mustKeepRanges)){
        if (length(intersect(allChromosomesMustKeep,part))>0){
          mustKeepWin <- getWinFromGranges(mustKeepRanges[seqnames(mustKeepRanges) %in% part],part,statInfo[idPart,],win,step)
        }
      }
      
      firstReadIndex <- ((floor(bam$flag/64) %% 2) == 1)
      secondReadIndex <- !firstReadIndex
      
      winFirstPositiveAlignments <- getWinOfAlignments(bam,"+",win,step,limit,firstReadIndex,coverage=coverage)
      winSecondPositiveAlignments <- getWinOfAlignments(bam,"+",win,step,limit,secondReadIndex,coverage=coverage)
      winFirstNegativeAlignments <- getWinOfAlignments(bam,"-",win,step,limit,firstReadIndex,coverage=coverage)
      winSecondNegativeAlignments <- getWinOfAlignments(bam,"-",win,step,limit,secondReadIndex,coverage=coverage)

      probaWinFirst <- keptProbaWin(winFirstPositiveAlignments,winFirstNegativeAlignments,win,step,logitThreshold,pvalueThreshold,errorRate,mustKeepWin,min,max,getWin,coverage=coverage) # the probability of each positive/negative window; this probability will be assigned to every positive/negative read in that window
      probaWinSecond <- keptProbaWin(winSecondPositiveAlignments,winSecondNegativeAlignments,win,step,logitThreshold,pvalueThreshold,errorRate,mustKeepWin,min,max,getWin,coverage=coverage) # the probability of each positive/negative window; this probability will be assigned to every positive/negative read in that window
      if (getWin){ 
        allWin <- rbind(allWin,data.frame("Type"="First",getWinInChromosome(probaWinFirst$Win,part,statInfo,win,step)),
                        data.frame("Type"="Second",getWinInChromosome(probaWinSecond$Win,part,statInfo,win,step)))
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
      for (i in seq_along(part)){
        id <- which(chromosomes == part[i])
        if (statInfo$NbOriginalReads[id]>0){
          rangeInChr <- which((keptAlignments>=statInfo$FirstReadInPartition[id] & keptAlignments<= statInfo$LastReadInPartition[id]))
          statInfo$NbOriginalFirstReads[id] <- sum(firstReadIndex[statInfo$FirstReadInPartition[id]:statInfo$LastReadInPartition[id]])
          statInfo$NbOriginalSecondReads[id] <- sum(secondReadIndex[statInfo$FirstReadInPartition[id]:statInfo$LastReadInPartition[id]])
          statInfo$NbKeptFirstReads[id] <- sum(firstReadIndex[rangeInChr])
          statInfo$NbKeptSecondReads[id] <- sum(secondReadIndex[rangeInChr])
        }
        if (statInfo$NbKeptFirstReads[id]>0 | statInfo$NbKeptSecondReads[i]>0){
          chromosomeIndex <- which(allChromosomes == part[i])
          bamSave(writer,bamRange(reader,c(chromosomeIndex-1,0,lengthSeq[chromosomeIndex]))[keptAlignments[rangeInChr]- statInfo$FirstReadInPartition[id]+1,],refid=chromosomeIndex-1)#write the kept alignments into the output bam file
          rm(rangeInChr)
        }
      }
      cat(paste0("Sequence ",part,", length: ",statInfo$Length[idPart],
                 ", number of first reads: ",statInfo$NbOriginalFirstReads[idPart],
                 ", number of second reads: ",statInfo$NbOriginalSecondReads[idPart],
                 ", number of kept first reads: ",statInfo$NbKeptFirstReads[idPart],
                 ", number of kept second reads: ",statInfo$NbKeptSecondReads[idPart],"\n"),file=statfile,append=append)
      if (!append) append <- TRUE
      rm(keptAlignments)
    }
  }
  bamClose(writer)
  bamClose(reader)
  cat("Summary:\n",file = statfile, append = append)
  cat(paste0("Number of original first reads: ",sum(statInfo$NbOriginalFirstReads),
             ", number of original second reads: ",sum(statInfo$NbOriginalSecondReads),
             ", number of kept first reads: ",sum(statInfo$NbKeptFirstReads),
             ", number of kept second reads: ",sum(statInfo$NbKeptSecondReads),
             ", removal proportion of first reads: ",(sum(statInfo$NbOriginalFirstReads)-sum(statInfo$NbKeptFirstReads))/sum(statInfo$NbOriginalFirstReads),
             ", removal proportion of second reads: ",(sum(statInfo$NbOriginalSecondReads)-sum(statInfo$NbKeptSecondReads))/sum(statInfo$NbOriginalSecondReads),"\n"),file = statfile,append=TRUE)
  endTime <- proc.time()
  cat(paste0("Total elapsed time ",(endTime-startTime)[[3]]/60," minutes\n"),file = statfile,append=TRUE)
  if (getWin){
    return(allWin %>% dplyr::mutate("Start" = (Start-1)*step+1))
  }
}



