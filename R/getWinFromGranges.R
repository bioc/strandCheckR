#' @title get the ranges of windows from a GRanges
#'
#' @export
#'
#'
getWinFromGranges <- function(mustKeepRanges,chromosomes,win,step){
  chroms <- intersect(seqnames(mustKeepRanges),chromosomes)
  if (length(chroms)>0){
    results <- lapply(chroms,function(chr){
      mkr <- mustKeepRanges[seqnames(mustKeepRanges)==chr]
      mustKeepPos <- getWinFromIRanges(ranges(mkr)[strand(mkr)!="-",],win,step,0) %>% coverage()
      mustKeepPos <- (mustKeepPos>0)
      mustKeepNeg <- getWinFromIRanges(ranges(mkr)[strand(mkr)!="+",],win,step,0) %>% coverage()
      mustKeepNeg <- (mustKeepNeg>0)
      return(list("Positive"=mustKeepPos,"Negative"=mustKeepNeg))
    })
    names(results) <- chroms
    return(results)
  }
  else{
    return(list())
  }
}
