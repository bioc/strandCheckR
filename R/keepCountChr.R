keepCountChr <- function(nbSample,alignments,chr,len,win,step,pvalueThreshold,minR,maxR,limit,logitThreshold){#compute the reads to be kept for the whole genome
  message("Chromosome: ",chr,", Length: ",len)
  position <- computePosition(alignments,chr)
  if (minR==0 && maxR==0 && !(missing(logitThreshold))){
    windows <- computeWinCount0(position$Pos$start,position$Pos$end,position$Neg$start,position$Neg$end,len,win,step,limit,logitThreshold) #compute information in each sliding windows
  }
  else{
    position <- reorder(position,win,step,limit)
    if (missing(logitThreshold)) 
      windows <- computeWinCountNoThreshold(position$Pos$start,position$Pos$end,position$Neg$start,position$Neg$end,len,win,step,minR,maxR,limit) #compute information in each sliding windows
    else 
      windows <- computeWinCount(position$Pos$start,position$Pos$end,position$Neg$start,position$Neg$end,len,win,step,minR,maxR,limit,logitThreshold) #compute information in each sliding windows
  }
  windows$Plus <- mutate(windows$Plus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% filter(pvalue<=pvalueThreshold) #get kept positive windows
  windows$Minus <- mutate(windows$Minus,"pvalue"=pnorm(value,lower.tail = FALSE)) %>% filter(pvalue<= pvalueThreshold) #get negative windows
  keptFrags <- keepRead(nbSample,position$Pos$start,position$Pos$end,position$Pos$group,position$Pos$sample,position$Neg$start,position$Neg$end,position$Neg$group,position$Neg$sample,windows$Plus$win,windows$Minus$win,len,win,step,limit) 
  remove(windows)
  gc()
  keep <- list()
  for (i in c(1:nbSample)){
    keep[[i]] <- c(position$Index[[i]]$Pos[unique(keptFrags$Pos[[i]])],position$Index[[i]]$Neg[unique(keptFrags$Neg[[i]])]) %>% sort() 
    message("Sample: ",i,", Number of reads: ",position$nbRead[i],", Number of kept reads: ",length(keep[[i]]))
  }
  rm(list=setdiff(ls(), "keep"))
  gc()
  return (keep)
}