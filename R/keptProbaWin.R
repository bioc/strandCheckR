#' @title get the probability of keeping each window
#'
#' @export
#' @importFrom IRanges coverage
#'
keptProbaWin <- function(winPositiveAlignments,winNegativeAlignments,win,step,logitThreshold,pvalueThreshold,errorRate,mustKeepWin,min,max,getWin,coverage=FALSE){
  if (coverage){
    lenPC <- length(winPositiveAlignments$Coverage)
    lenNC <- length(winNegativeAlignments$Coverage)
    if (lenPC>lenNC) {winNegativeAlignments$Coverage <- c(winNegativeAlignments$Coverage,rep(0,lenPC-lenNC))} else {winPositiveAlignments$Coverage <- c(winPositiveAlignments$Coverage,rep(0,lenNC-lenPC))}
    nbWin <- ceiling((length(winPositiveAlignments$Coverage)-win)/step)+1
    nbPositiveReads <- Views(winPositiveAlignments$Coverage,
                             start = seq(1,(nbWin-1)*step+1,step),
                             end=seq(win,(nbWin-1)*step+win,step)) %>%
      sum() %>% Rle()
    nbNegativeReads <- Views(winNegativeAlignments$Coverage,
                             start = seq(1,(nbWin-1)*step+1,step),
                             end=seq(win,(nbWin-1)*step+win,step)) %>%
      sum() %>% Rle()
    if (min>0 | max>0 | getWin){
      maxCoverage <- Views(winPositiveAlignments$Coverage+winNegativeAlignments$Coverage,
                           start = seq(1,(nbWin-1)*step+1,step),
                           end=seq(win,(nbWin-1)*step+win,step)) %>% 
        max() %>% Rle()
    }
  }
  else{
    nbPositiveReads <- coverage(winPositiveAlignments$Win)
    nbNegativeReads <- coverage(winNegativeAlignments$Win)
    lenP <- length(nbPositiveReads)
    lenN <- length(nbNegativeReads)
    if (lenP>lenN) { nbNegativeReads <- c(nbNegativeReads,rep(0,lenP-lenN))} else {nbPositiveReads <- c(nbPositiveReads,rep(0,lenN-lenP))}
  }
  toTest <- (logitThreshold-abs(log(nbPositiveReads/nbNegativeReads)))/sqrt((nbPositiveReads+nbNegativeReads)/nbPositiveReads/nbNegativeReads)
  pvalue <- Rle(pnorm(runValue(toTest)),runLength(toTest))
  rm(toTest)
  pvalue[(nbPositiveReads==0 | nbNegativeReads==0)] <- 0
  if (min>0){
    if (coverage){
      pvalue[maxCoverage<min] <- 1
    }
    else{
      pvalue[(nbPositiveReads+nbNegativeReads)<min] <- 1
    }
  }
  keptWin <- rep(TRUE,length(runValue(pvalue)))
  keptWin[runValue(pvalue)>pvalueThreshold] <- FALSE
  keptWin <- Rle(keptWin,runLength(pvalue))

  keptProbaPosWin <- keptWin*((nbPositiveReads>nbNegativeReads)*((nbPositiveReads-nbNegativeReads)/nbPositiveReads)+(nbPositiveReads<nbNegativeReads)*errorRate)
  keptProbaNegWin <- keptWin*((nbPositiveReads<nbNegativeReads)*((nbNegativeReads-nbPositiveReads)/nbNegativeReads)+(nbPositiveReads>nbNegativeReads)*errorRate)
  
  if (length(mustKeepWin)>0){
    if (length(mustKeepWin$Positive)>length(keptProbaPosWin)) mustKeepWin$Positive <- mustKeepWin$Positive[length(keptProbaPosWin)]
    else mustKeepWin$Positive <- c(mustKeepWin$Positive,rep(0,length(keptProbaPosWin)-length(mustKeepWin$Positive)))
    keptProbaPosWin <- mustKeepWin$Positive + (!mustKeepWin$Positive)*keptProbaPosWin
    if (length(mustKeepWin$Negative)>length(keptProbaNegWin)) mustKeepWin$Negative <- mustKeepWin$Negative[length(keptProbaNegWin)]
    else mustKeepWin$Negative <- c(mustKeepWin$Negative,rep(0,length(keptProbaNegWin)-length(mustKeepWin$Negative)))
    keptProbaNegWin <- mustKeepWin$Negative + (!mustKeepWin$Negative)*keptProbaNegWin
  }
  if (max>0){
    if (coverage){
      keepMorePos <- (nbPositiveReads>nbNegativeReads)*(maxCoverage>=max)
      keepMoreNeg <- (nbPositiveReads<nbNegativeReads)*(maxCoverage>=max)
    }
    else{
      keepMorePos <- (nbPositiveReads>nbNegativeReads)*(nbPositiveReads>=max)
      keepMoreNeg <- (nbPositiveReads<nbNegativeReads)*(nbNegativeReads>=max)
    }
    keptProbaPosWin <- keepMorePos+(!keepMorePos)*keptProbaPosWin
    keptProbaNegWin <- keepMoreNeg+(!keepMoreNeg)*keptProbaNegWin
  }
  if (getWin){
      presentWin <- which(as.vector((nbPositiveReads>0) | (nbNegativeReads>0)) == TRUE)
      if (length(presentWin)>0){
        if (coverage){
          Win <- (data.frame("Start" = presentWin, "NbPositiveReads" = nbPositiveReads[presentWin], "NbNegativeReads" = nbNegativeReads[presentWin],"MaxCoverage" = maxCoverage[presentWin])) 
        }
        else{
          Win <- (data.frame("Start" = presentWin, "NbPositiveReads" = nbPositiveReads[presentWin], "NbNegativeReads" = nbNegativeReads[presentWin]))  
        }
      }
    return(list("Positive"=keptProbaPosWin,"Negative"=keptProbaNegWin,"Win"=Win))
  }
  else{
    return(list("Positive"=keptProbaPosWin,"Negative"=keptProbaNegWin))
  }
}


