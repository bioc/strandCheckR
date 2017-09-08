#' @title get the probability of begin kept for each window
#' @description calculate the keeping probability of each window based on its positive/negative proportion
#' @param winPositiveAlignments an object returned by getWinOfAlignments for positive reads
#' @param winNegativeAlignments an object returned by getWinOfAlignments for negative reads
#' @param winWidth the width of the sliding window, 1000 by default.
#' @param winStep the winStep length to sliding the window, 100 by default.
#' @param logitThreshold logistic value of the threshold
#' @param pvalueThreshold threshold of p-value
#' @param min In the case that \code{useCoverage=FALSE}, if a window has least than \code{min} reads, then it will be rejected regardless the strand proportion. 
#'        For the case that \code{useCoverage=TRUE}, if a window has max coverage least than \code{min}, then it will be rejected. 0 by default
#' @param max In the case that \code{useCoverage=FALSE},if a window has more than \code{max} reads, then it will be kept regardless the strand proportion. 
#'        For the case that \code{useCoverage=TRUE}, if a window has max coverage more than \code{max}, then it will be kept. 
#'        If 0 then it doesn't have effect on selecting window. 0 by default.
#' @param errorRate the probability that an RNA read takes the false strand. 0.01 by default
#' @param useCoverage if TRUE, then the strand information in each window corresponds to the sum of coverage coming from positive/negative reads; and not the number of positive/negative reads as default.
#' @param mustKeepWin the windows that must be kept regardless their strand proportion
#' @param getWin if TRUE, the function will return a data frame containing the information of all windows. It's FALSE by default.
#' 
#' @importFrom IRanges coverage Views
#' @importFrom stats pnorm
#'
#' @export
keptProbaWin <- function(winPositiveAlignments,winNegativeAlignments,winWidth,winStep,logitThreshold,pvalueThreshold,errorRate,mustKeepWin,min,max,getWin,useCoverage=FALSE){
  if (useCoverage){
    fromCoverage <- calculateStrandCoverage(winPositiveAlignments,winNegativeAlignments,winWidth,winStep)
    nbPositive <- fromCoverage$CovPositive
    nbNegative <- fromCoverage$CovNegative
  }
  else{
    fromNbReads <- calculateStrandNbReads(winPositiveAlignments,winNegativeAlignments)
    nbPositive <- fromNbReads$NbPositive
    nbNegative <- fromNbReads$NbNegative
  }
  
  if (getWin){
    presentWin <- which(as.vector((nbPositive>0) | (nbNegative>0)) == TRUE)
    if (length(presentWin)>0){
      Win <- DataFrame(Start = presentWin, 
                       CovPositive = nbPositive[presentWin], CovNegative = nbNegative[presentWin])
      if (useCoverage){
        Win[["MaxCoverage"]] = fromCoverage$maxCoverage[presentWin]
      }
    }
  }
  
  toTest <- (logitThreshold-abs(log(nbPositive/nbNegative)))/sqrt((nbPositive+nbNegative)/nbPositive/nbNegative)
  pvalue <- Rle(pnorm(runValue(toTest)),runLength(toTest))
  rm(toTest)
  pvalue[(nbPositive==0 | nbNegative==0)] <- 0
  if (min>0){
    if (useCoverage){
      pvalue[fromCoverage$MaxCoverage<min] <- 1
    }
    else{
      pvalue[(nbPositive+nbNegative)<min] <- 1
    }
  }
  keptWin <- rep(TRUE,length(runValue(pvalue)))
  keptWin[runValue(pvalue)>pvalueThreshold] <- FALSE
  keptWin <- Rle(keptWin,runLength(pvalue))

  keptProbaPosWin <- keptWin*((nbPositive>nbNegative)*((nbPositive-nbNegative)/nbPositive)+(nbPositive<nbNegative)*errorRate)
  keptProbaNegWin <- keptWin*((nbPositive<nbNegative)*((nbNegative-nbPositive)/nbNegative)+(nbPositive>nbNegative)*errorRate)
  
  if (length(mustKeepWin)>0){
    if (length(mustKeepWin$Positive)>length(keptProbaPosWin)) mustKeepWin$Positive <- mustKeepWin$Positive[length(keptProbaPosWin)]
    else mustKeepWin$Positive <- c(mustKeepWin$Positive,rep(0,length(keptProbaPosWin)-length(mustKeepWin$Positive)))
    keptProbaPosWin <- mustKeepWin$Positive + (!mustKeepWin$Positive)*keptProbaPosWin
    if (length(mustKeepWin$Negative)>length(keptProbaNegWin)) mustKeepWin$Negative <- mustKeepWin$Negative[length(keptProbaNegWin)]
    else mustKeepWin$Negative <- c(mustKeepWin$Negative,rep(0,length(keptProbaNegWin)-length(mustKeepWin$Negative)))
    keptProbaNegWin <- mustKeepWin$Negative + (!mustKeepWin$Negative)*keptProbaNegWin
  }
  if (max>0){
    if (useCoverage){
      keepMorePos <- (nbPositive>nbNegative)*(fromCoverage$MaxCoverage>=max)
      keepMoreNeg <- (nbPositive<nbNegative)*(fromCoverage$MaxCoverage>=max)
    }
    else{
      keepMorePos <- (nbPositive>nbNegative)*(nbPositive>=max)
      keepMoreNeg <- (nbPositive<nbNegative)*(nbNegative>=max)
    }
    keptProbaPosWin <- keepMorePos+(!keepMorePos)*keptProbaPosWin
    keptProbaNegWin <- keepMoreNeg+(!keepMoreNeg)*keptProbaNegWin
  }
  if (getWin){
    return(list("Positive"=keptProbaPosWin,"Negative"=keptProbaNegWin,"Win"=Win))
  }
  else{
    return(list("Positive"=keptProbaPosWin,"Negative"=keptProbaNegWin))
  }
}


