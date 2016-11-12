#' @title Calculate the Reads to be Kept for one Chromosome
#' 
#' @description None
#' 
#' @param nbSample number of samples in consideration
#' @param alignmentChr the alignments of the considerint chromosome for all samples
#' @param len the length of the considering chromosome
#' @param win the length of the sliding window
#' @param step the step length to slide the window
#' @param pvalueThreshold the threshold for the p-value
#' @param minCov if a window has the max coverage least than minCov, then it will be rejected
#' @param maxCov if a window has the max coverage greater than maxCov, then it will be kept
#' @param limit the proportion that a read must be inside a window in order to be counted
#' @param threshold the threshold upper which we keep the reads
#' 

keepCountChr <- function(nbSample,alignmentChr,len,win,step,pvalueThreshold,minCov,maxCov,limit,threshold){#compute the reads to be kept for the whole genome
  logitThreshold <- binomial()$linkfun(threshold)
  position <- computePosition(alignmentChr)
  if (minCov==0 && maxCov==0){
    windows <- computeWinCount0(position$Pos$start,position$Pos$end,position$Neg$start,position$Neg$end,len,win,step,limit,logitThreshold) #compute information in each sliding windows
  }
  else{
    position <- reorder(position,win,step,limit)
    windows <- computeWinCount(position$Pos$start,position$Pos$end,position$Neg$start,position$Neg$end,len,win,step,minCov,maxCov,limit,logitThreshold) #compute information in each sliding windows
  }
  windows$Plus <- mutate(windows$Plus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% filter(pvalue<=pvalueThreshold) #get kept positive windows
  windows$Minus <- mutate(windows$Minus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% filter(pvalue<= pvalueThreshold) #get kept negative windows
  keptFrags <- keepRead(length(alignmentChr),position$Pos$start,position$Pos$end,position$Pos$group,position$Pos$sample,position$Neg$start,position$Neg$end,position$Neg$group,position$Neg$sample,windows$Plus$win,windows$Minus$win,len,win,step,limit) 
  remove(windows)
  keep <- list()
  for (i in seq_along(alignmentChr)){
    keep[[i]] <- c(position$Index[[i]]$Pos[unique(keptFrags$Pos[[i]])],position$Index[[i]]$Neg[unique(keptFrags$Neg[[i]])]) %>% sort() 
    message("Sample: ",i,", Number of reads: ",position$nbRead[i],", Number of kept reads: ",length(keep[[i]]))
  }
  rm(list=setdiff(ls(), "keep"))
  gc()
  return (keep)
}