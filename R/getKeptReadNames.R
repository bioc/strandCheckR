
getKeptReadNames <- function(alignment,covPos,covNeg,mustKeepPos,mustKeepNeg,getWin,readLength,len,win,step,pvalueThreshold,minCov,maxCov,logitThresholdP,logitThresholdM,errorRate){
  windows <- computeWinInfo(runLength(covPos),runValue(covPos),runLength(covNeg),runValue(covNeg),readLength,len,win,step,minCov)
  
  plus <- dplyr::filter(windows,NbPositiveReads>=NbNegativeReads) 
  if (nrow(plus)>0){
    plus <- dplyr::mutate(plus,"propor" = NbPositiveReads/(NbPositiveReads+NbNegativeReads),"nbReads"=NbPositiveReads+NbNegativeReads)
    pvalueP <- pnorm(logitThresholdP,mean=binomial()$linkfun(plus$propor),sd=sqrt(1/(plus$nbReads)/plus$propor/(1-plus$propor)))
    plus <- dplyr::mutate(plus,"pvalue"=pvalueP) %>% 
      dplyr::filter(pvalueP<=pvalueThreshold | Start %in% mustKeepPos | MaxCoverage >= maxCov) 
    if (nrow(plus)>0){
      plus <- dplyr::mutate(plus,"Start" = floor(Start/step)+1) %>% 
        dplyr::select(c(Start,propor))
    }
  }
  minus <- dplyr::filter(windows,NbPositiveReads<NbNegativeReads) 
  if (nrow(minus)>0){
    minus <- dplyr::mutate(minus,"propor" = NbPositiveReads/(NbPositiveReads+NbNegativeReads),"nbReads"=NbPositiveReads+NbNegativeReads)
    pvalueM <- pnorm(logitThresholdM,mean=binomial()$linkfun(minus$propor),sd=sqrt(1/(minus$nbReads)/minus$propor/(1-minus$propor)),lower.tail = FALSE)
    minus <- dplyr::mutate(minus,"pvalue"=pvalueM) %>% 
      dplyr::filter(pvalueM<=pvalueThreshold | Start %in% mustKeepNeg | MaxCoverage >= maxCov)
    if (nrow(minus)>0){
      minus <- dplyr::mutate(minus,"Start" = floor(Start/step)+1) %>% 
        dplyr::select(c(Start,propor)) 
    }
  } 
  nameReads <- c()  
  if (nrow(plus)>0 || nrow(minus)>0){
    strand <- as.vector(strand(alignment))
    indexPos <- which(strand=="+")
    indexNeg <- which(strand=="-")
    remove(strand)
    fragments <- getFragment(alignment)
    keptReads <- keepRead(fragments$Pos,fragments$Neg,plus,minus,win,step,errorRate);   
    keptReads <- c(indexPos[unique(keptReads$Pos)],indexNeg[unique(keptReads$Neg)]) %>% sort() 
    nameReads <- names(alignment)[keptReads] %>% unique()
  }
  if (getWin){
    return(list("Win"=windows,"nameReads"=nameReads))
  }
  else{return(nameReads)}
}
