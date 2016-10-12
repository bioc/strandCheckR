#' filter several bamfiles based on the strand specific of all of them. Strand of sliding windows are caculated based on read counts.
#' @param bamfilein the input bam files to be filterd
#' @param bamfileout the output bam files
#' @param chromosomes the list of chromosomes to be filtered
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' @param threshold the threshold upper which we keep the reads
#' @param pvalueThreshold the threshold for the p-value
#' @export
filterMultiCount <- function(bamfilein,bamfileout,chromosomes,win,step,threshold,pvalueThreshold,minR){
  library(GenomicAlignments)
  library(rbamtools)
  library(Rcpp)
  library(dplyr)
  logitThreshold <- binomial()$linkfun(threshold)
  alignments <- list()
  reader <- list()
  header <- list()
  writer <- list()
  for (i in c(1:length(bamfilein))){
    alignments[[i]] <- readGAlignments(bamfilein[i],param=ScanBamParam(what=c("cigar")))#read file i
    reader[[i]] <- bamReader(bamfilein[i],idx=TRUE) #open a reader of the input bamfile to extract read afterward
    header[[i]] <- getHeader(reader[[i]]) #get the header of the input bam file
    writer[[i]] <- bamWriter(header[[i]],bamfileout[i]) #prepare to write the output bamfile with the same header
  }
  refSeqs <- getRefData(reader[[1]])
  allChromosomes <- refSeqs$SN
  lenSeq <- refSeqs$LN
  if (is.null(chromosomes)) chromosomes<-allChromosomes
  for (chr in chromosomes){
    message("Chromosome ",chr)
    positionPos <- data.frame("group"=c(),"start"=c(),"end"=c(),"sample"=c()) #position of each fragment of positive reads
    positionNeg <- data.frame("group"=c(),"start"=c(),"end"=c(),"sample"=c()) #position of each fragment of negative reads
    index <- list() #the index of each read in the alignment from its index in positive/negative read lists
    nbRead <- c() #number of original reads 
    for (i in c(1:length(bamfilein))){ #compute positionPos and positionNeg for all samples
      al <- alignments[[i]][seqnames(alignments[[i]])==chr]
      nbRead <- c(nbRead,length(al))
      index[[i]] <- getIndex(as.vector(strand(al)))
      alPos<-al[strand(al)=="+"]
      pos <- extractAlignmentRangesOnReference(cigar(alPos),pos=start(alPos)) %>% data.frame() %>% dplyr::select(-c(group_name,width))
      remove(alPos)
      pos[["sample"]] <- rep(i,nrow(pos))
      positionPos <- rbind(positionPos,pos)
      remove(pos)
      alNeg<-al[strand(al)=="-"]
      remove(al)
      neg <- extractAlignmentRangesOnReference(cigar(alNeg),pos=start(alNeg)) %>% data.frame() %>% dplyr::select(-c(group_name,width))
      remove(alNeg)
      neg[["sample"]] <- rep(i,nrow(neg))
      positionNeg <- rbind(positionNeg,neg)
      remove(neg)
    } 
    positionPos[["firstW"]] <- ceiling((positionPos$start-win-1+floor((positionPos$end-positionPos$start+1)*25/100))/step)
    positionNeg[["firstW"]] <- ceiling((positionNeg$start-win-1+floor((positionNeg$end-positionNeg$start+1)*25/100))/step)
    positionPos <-  positionPos[order(positionPos$firstW),] #reorder positionPos following the starting position of each fragment
    positionNeg <-  positionNeg[order(positionNeg$firstW),] #reorder positionNeg following the starting position of each fragment
    chromosomeIndex <- which(chromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    windows <- computeWinCount(positionPos$start,positionPos$end,positionNeg$start,positionNeg$end,len,win,step,logitThreshold,minR) #compute information in each sliding windows
    windows$Plus["pvalue"] <- pnorm(windows$Plus$value,lower.tail = FALSE) #compute pvalue for positive windows
    windows$Minus["pvalue"] <- pnorm(windows$Minus$value,lower.tail = FALSE)#compute pvalue for negative windows
    keepWinPos <- filter(windows$Plus, pvalue <= pvalueThreshold)$win # the indices of positive windows to be kept
    keepWinNeg <- filter(windows$Minus, pvalue <= pvalueThreshold)$win # the indices of negative windows to be kept
    remove(windows)
    reads <- keepRead(length(bamfilein),positionPos$start,positionPos$end,positionPos$group,positionPos$sample,positionNeg$start,positionNeg$end,positionNeg$group,positionNeg$sample,keepWinPos,keepWinNeg,len,win,step) #compute index of positive/negative reads to be kept
    remove(positionPos)
    remove(positionNeg)
    remove(keepWinPos)
    remove(keepWinNeg)
    for (i in c(1:length(bamfilein))){
      keep <- c(index[[i]]$Pos[unique(reads$Pos[[i]])],index[[i]]$Neg[unique(reads$Neg[[i]])]) %>% sort() #index of reads to be kept
      message("Sample ",i,", number of reads: ",nbRead[i],", number of kept reads: ",length(keep))
      if (length(keep)>0){
        #get the range of kept reads
        range <- bamRange(reader[[i]],c(chromosomeIndex-1,0,lenSeq[chromosomeIndex]))
        #write the kept reads into output file
        bamSave(writer[[i]],range[keep,],refid=chromosomeIndex-1)
        remove(range)
      }
      remove(keep)
    }
    remove(index)
    remove(reads)
  }
  remove(alignments)
  for (i in (1:length(bamfilein))) bamClose(writer[[i]])
}

