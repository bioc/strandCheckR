#' @title get the ranges of windows from a GRanges
#' @description get the lists of positive/negative windows that overlap a GRanges objects
#' @param gr a GRanges object
#' @param chromosomes a list of chromosome names
#' @param statinfo a data frame that contains some information of the alignments
#' @param win the length of sliding window
#' @param step the step length to sliding the window
#' @export
#'
#'
getWinFromGranges <- function(gr,chromosomes,statinfo,win,step){
  for (i in seq_along(chromosomes)){
    r <- which(as.vector(seqnames(gr)) == chromosomes[i])
    if (length(r)>0){
      start(ranges(gr)[r]) <- start(ranges(gr)[r]) + statinfo$FirstBaseInPartition[i] -1
      end(ranges(gr)[r]) <- end(ranges(gr)[r]) + statinfo$FirstBaseInPartition[i] -1  
    }
  }
  mustKeepPos <- getWinFromIRanges(ranges(gr)[strand(gr)!="-",],win,step,1) %>% coverage()
  mustKeepPos <- (mustKeepPos>0)
  mustKeepNeg <- getWinFromIRanges(ranges(gr)[strand(gr)!="+",],win,step,1) %>% coverage()
  mustKeepNeg <- (mustKeepNeg>0)
  return(list("Positive"=mustKeepPos,"Negative"=mustKeepNeg))
}

