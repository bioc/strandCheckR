
getKeptReadNames <- function(alignment,covPos,covNeg,mustKeepPosWin,mustKeepNegWin,getWin,readLength,len,win,step,pvalueThreshold,minCov,maxCov,logitThresholdP,logitThresholdM,errorRate){
  if (getWin){#get details of each window for filtering and plotting
    windows <- computeWinVerbose(runLength(covPos),runValue(covPos),runLength(covNeg),runValue(covNeg),readLength,len,win,step,minCov,maxCov)
  }
  else{#get details of each window for filtering
    windows <- computeWin(runLength(covPos),runValue(covPos),runLength(covNeg),runValue(covNeg),readLength,len,win,step,minCov,maxCov)
  }
  
  windows$Plus$pvalue <- pnorm(logitThresholdP,mean=binomial()$linkfun(windows$Plus$propor),sd=windows$Plus$error)
  windows$Plus <- dplyr::filter(windows$Plus,pvalue<=pvalueThreshold | win %in% mustKeepPosWin) %>%
    dplyr::select(c(win,propor)) 
  windows$Minus$pvalue <- pnorm(logitThresholdM,mean=binomial()$linkfun(windows$Minus$propor),sd=windows$Minus$error,lower.tail = FALSE)
  windows$Minus <- dplyr::filter(windows$Minus,pvalue<=pvalueThreshold | win %in% mustKeepNegWin) %>% 
    dplyr::select(c(win,propor)) 
  
  nameReads <- c()  
  if (nrow(windows$Plus)>0 || nrow(windows$Minus)>0){
    strand <- as.vector(strand(alignment))
    indexPos <- which(strand=="+")
    indexNeg <- which(strand=="-")
    remove(strand)
    fragments <- getFragment(alignment)
    keptReads <- keepRead(fragments$Pos,fragments$Neg,windows$Plus,windows$Minus,win,step,errorRate);   
    keptReads <- c(indexPos[unique(keptReads$Pos)],indexNeg[unique(keptReads$Neg)]) %>% sort() 
    nameReads <- names(alignment)[keptReads] %>% unique()
  }
  if (getWin){
    return(list("Win"=windows$Win,"nameReads"=nameReads))
  }
  else{return(nameReads)}
}
