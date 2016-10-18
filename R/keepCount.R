keepCount <- function(bamfilein,bamfileout,chromosomes,allChromosomes,lenSeq,win,step,threshold,pvalueThreshold,minR,maxR,limit){#compute the reads to be kept for the whole genome
  logitThreshold <- binomial()$linkfun(threshold)
  alignments <- list()
  for (i in c(1:length(bamfilein))){
    alignments[[i]] <- readGAlignments(bamfilein[i],param=ScanBamParam(what=c("cigar")))#read file i
  }
  keepReads <- list()
  for (f in c(1:length(bamfilein))){
    keepReads[[i]] <- list()
  }
  for (chr in chromosomes){
    chromosomeIndex <- which(allChromosomes==chr)
    len <- lenSeq[chromosomeIndex]
    message("Chromosome: ",chr,", Length: ",len)
    position <- computePosition(alignments,chr)
    if (minR==0 && maxR==0){
      windows <- computeWinCount0(position$Pos$start,position$Pos$end,position$Neg$start,position$Neg$end,len,win,step,logitThreshold,maxR,limit) #compute information in each sliding windows
    }
    else{
      position <- reorder(position,win,step,limit)
      windows <- computeWinCount(position$Pos$start,position$Pos$end,position$Neg$start,position$Neg$end,len,win,step,logitThreshold,minR,maxR,limit) #compute information in each sliding windows
    }
    windows$Plus <- mutate(windows$Plus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% filter(pvalue<=pvalueThreshold) #compute pvalue for positive windows
    windows$Minus <- mutate(windows$Minus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% filter(pvalue<= pvalueThreshold) #compute pvalue for negative windows
    #message("keep win Pos: ",nrow(windows$Plus))
    #message("keep win Neg: ",nrow(windows$Minus))
    #compute index of positive/negative reads to be kept
    keptFrags <- keepRead(length(bamfilein),position$Pos$start,position$Pos$end,position$Pos$group,position$Pos$sample,position$Neg$start,position$Neg$end,position$Neg$group,position$Neg$sample,windows$Plus$win,windows$Minus$win,len,win,step,limit) 
    remove(windows)
    gc()
    for (i in c(1:length(bamfilein))){
      #message("keep frag Pos: ",length(keptFrags$Pos[[i]]))
      #message("keep frag Neg: ",length(keptFrags$Neg[[i]]))
      keepReads[[i]][[chr]] <- c(position$Index[[i]]$Pos[unique(keptFrags$Pos[[i]])],position$Index[[i]]$Neg[unique(keptFrags$Neg[[i]])]) %>% sort() #index of reads to be kept
      message("Sample: ",i,", Number of reads: ",position$nbRead[i],", Number of kept reads: ",length(keepReads[[i]][[chr]]))
    }
    remove(keptFrags)
    gc()
    remove(position)
    gc()
  }
  remove(alignments)
  gc()
  return (keepReads)
}