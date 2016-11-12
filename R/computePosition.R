#' @title Calculate Alignment Positions
#' 
#' @description Calculate start, end position and strand of every alignment of the considering chromosome. 
#' Each alignment is a group, consisting of all of its splited fragments (if an alignment is not splitted then its group contain only one fragment).
#' @param bamfilein the input bam files to be filterd
#' @param chromosomes the list of chromosomes to be filtered
#' @param allChromosomes the list of all chromosomes
#' @param lenSeq the length of every chromosome
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' @param pvalueThreshold the threshold for the p-value
#' @param minCov if a window has the max coverage least than minCov, then it will be rejected
#' @param maxCov if a window has the max coverage greater than maxCov, then it will be kept
#' @param limit the proportion that a read must be inside a window in order to be counted
#' @param threshold the threshold upper which we keep the reads
#' 
#' @return the information about the reads
#' nbRead: total number of read
#' positionPos: data frame of fragments come from positive reads
#' positionNeg: data frame of fragments come from negative reads
#' index: the indices of position/negative reads in the whole alignments

computePosition <- function(alignmentChr){
  positionPos <- data.frame("group"=c(),"start"=c(),"end"=c(),"sample"=c()) #position of each fragment of positive reads
  positionNeg <- data.frame("group"=c(),"start"=c(),"end"=c(),"sample"=c()) #position of each fragment of negative reads
  index <- list() #the index of each read in the alignment from its index in positive/negative read lists
  nbRead <- c() #the number of reads
  for (i in c(1:length(alignmentChr))){ #compute positionPos and positionNeg for all samples
    nbRead <- c(nbRead,length(alignmentChr[[i]]))
    index[[i]] <- getIndex(as.vector(strand(alignmentChr[[i]])))
    alPos<-alignmentChr[[i]][strand(alignmentChr[[i]])=="+"]
    positionPos <- rbind(positionPos,extractAlignmentRangesOnReference(cigar(alPos),pos=start(alPos)) 
                         %>% data.frame() 
                         %>% dplyr::select(-c(group_name,width)) 
                         %>% mutate(sample=i))
    remove(alPos)
    alNeg<-alignmentChr[[i]][strand(alignmentChr[[i]])=="-"]
    positionNeg <- rbind(positionNeg,extractAlignmentRangesOnReference(cigar(alNeg),pos=start(alNeg)) 
                         %>% data.frame() 
                         %>% dplyr::select(-c(group_name,width))
                         %>% mutate(sample=i))
    remove(alNeg)
  } 
  rm(list=setdiff(ls(), c("nbRead","positionPos","positionNeg","index")))
  gc()
  return (list(nbRead=nbRead,Pos=positionPos,Neg=positionNeg,Index=index))
}
