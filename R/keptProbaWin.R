#' @title get the probability of begin kept for each window
#' @description calculate the keeping probability of each window based on its 
#' positive/negative proportion
#' @param winPositiveAlignments an object returned by getWinOfAlignments for 
#' positive reads
#' @param winNegativeAlignments an object returned by getWinOfAlignments for 
#' negative reads
#' @param winWidth the width of the sliding window, 1000 by default.
#' @param winStep the winStep length to sliding the window, 100 by default.
#' @param logitThreshold logistic value of the threshold
#' @param pvalueThreshold threshold of p-value
#' @param minCov In the case that \code{useCoverage=FALSE}, if a window has 
#' least than \code{minCov} reads, then it will be rejected regardless the 
#' strand proportion. 
#' For the case that \code{useCoverage=TRUE}, if a window has max coverage
#' least than \code{minCov}, then it will be rejected. 0 by default
#' @param maxCov In the case that \code{useCoverage=FALSE}, if a window has 
#' more than \code{maxCov} reads, then it will be kept regardless the strand 
#' proportion. 
#' For the case that \code{useCoverage=TRUE}, if a window has max coverage 
#' more than \code{maxCov}, then it will be kept. If 0 then it doesn't have 
#' effect on selecting window. 0 by default.
#' @param errorRate the probability that an RNA read takes the false strand. 
#' 0.01 by default
#' @param useCoverage if TRUE, then the strand information in each window 
#' corresponds to the sum of coverage coming from positive/negative reads; 
#' and not the number of positive/negative reads as default.
#' @param mustKeepWin the windows that must be kept regardless their strand 
#' proportion
#' @param getWin if TRUE, the function will return a data frame containing the 
#' information of all windows. It's FALSE by default.
#' 
#' @return A list of 2 numeric-Rle objects containing keeping probability of 
#' each +/- alignments. 
#' If \code{getWin=TRUE} then the list contains an additional DataFrame for the 
#' number of reads and coverage of the input window +/- alignments
#' 
#' @importFrom IRanges coverage Views
#' @importFrom stats pnorm
#'

keptProbaWin <- function(winPositiveAlignments, winNegativeAlignments, 
                        winWidth, winStep, logitThreshold, pvalueThreshold, 
                        errorRate, mustKeepWin, minCov, maxCov, getWin, 
                        useCoverage=FALSE){
    if (getWin){
        fromCoverage <- calculateStrandCoverage(winPositiveAlignments, 
                                                winNegativeAlignments, 
                                                winWidth,winStep)
        fromNbReads <- calculateStrandNbReads(winPositiveAlignments,
                                                winNegativeAlignments)
        presentWin <- which(as.vector(fromCoverage$CovPositive>0 | 
                                        fromCoverage$CovNegative>0) == TRUE)
        Win <- DataFrame(Type = "", Seq = "", 
                        Start = presentWin, End = 0,
                        NbPositive = fromNbReads$NbPositive[presentWin], 
                        NbNegative = fromNbReads$NbNegative[presentWin],
                        CovPositive = fromCoverage$CovPositive[presentWin], 
                        CovNegative = fromCoverage$CovNegative[presentWin],
                        MaxCoverage = fromCoverage$MaxCoverage[presentWin],
                        File = "")
    } else if (useCoverage){
        fromCoverage <- calculateStrandCoverage(winPositiveAlignments, 
                                                winNegativeAlignments, 
                                                winWidth,winStep)
    } else{
        fromNbReads <- calculateStrandNbReads(winPositiveAlignments,
                                            winNegativeAlignments)
    }
    if (useCoverage){
        pos <- fromCoverage$CovPositive
        neg <- fromCoverage$CovNegative
    } else{
        pos <- fromNbReads$NbPositive
        neg <- fromNbReads$NbNegative
    }
    toTest <- (logitThreshold-abs(log(pos/neg)))/sqrt((pos+neg)/pos/neg)
    pvalue <- Rle(pnorm(runValue(toTest)),runLength(toTest))
    rm(toTest)
    pvalue[(pos==0 | neg==0)] <- 0
    if (minCov>0){
        if (useCoverage){
            pvalue[fromCoverage$MaxCoverage<minCov] <- 1
        } else{
            pvalue[(pos+neg)<minCov] <- 1
        }
    }
    keptWin <- rep(TRUE,length(runValue(pvalue)))
    keptWin[runValue(pvalue)>pvalueThreshold] <- FALSE
    keptWin <- Rle(keptWin,runLength(pvalue))

    keptProbaPosWin <- keptWin*((pos>neg)*((pos-neg)/pos)+(pos<neg)*errorRate)
    keptProbaNegWin <- keptWin*((pos<neg)*((neg-pos)/neg)+(pos>neg)*errorRate)

    if (length(mustKeepWin)>0){
        if (length(mustKeepWin$Positive)>length(keptProbaPosWin)){ 
            mustKeepWin$Positive <- mustKeepWin$
                                            Positive[length(keptProbaPosWin)]
        } else {
            mustKeepWin$Positive <- c(mustKeepWin$Positive,
                    rep(0,length(keptProbaPosWin)-length(mustKeepWin$Positive)))
        }
        keptProbaPosWin <- mustKeepWin$Positive + 
            (!mustKeepWin$Positive)*keptProbaPosWin
        if (length(mustKeepWin$Negative)>length(keptProbaNegWin)){ 
            mustKeepWin$Negative <- mustKeepWin$
                                            Negative[length(keptProbaNegWin)]
        } else {
            mustKeepWin$Negative <- c(mustKeepWin$Negative,
                    rep(0,length(keptProbaNegWin)-length(mustKeepWin$Negative)))
        }
        keptProbaNegWin <- mustKeepWin$Negative + 
            (!mustKeepWin$Negative)*keptProbaNegWin
    }
    if (maxCov>0){
        if (useCoverage){
            keepMorePos <- (pos>neg)*(fromCoverage$MaxCoverage>=maxCov)
            keepMoreNeg <- (pos<neg)*(fromCoverage$MaxCoverage>=maxCov)
        }
        else{
            keepMorePos <- (pos>neg)*(pos>=maxCov)
            keepMoreNeg <- (pos<neg)*(neg>=maxCov)
        }
        keptProbaPosWin <- keepMorePos+(!keepMorePos)*keptProbaPosWin
        keptProbaNegWin <- keepMoreNeg+(!keepMoreNeg)*keptProbaNegWin
    }
    if (getWin){
        return(list("Positive"=keptProbaPosWin,
                    "Negative"=keptProbaNegWin,"Win"=Win))
    } else{
        return(list("Positive"=keptProbaPosWin,"Negative"=keptProbaNegWin))
    }
}


