#' @title get the ranges of windows from on IRanges
#' @description get the list of windows that overlaps an Iranges object
#' @param gr an Iranges object
#' @param win the length of sliding window
#' @param step the step length to sliding the window
#' @param limit a read is considered to be included in a window if and only if at least limit percent of it is in the window. 
#' 
#' @export
#'
#'
getWinFromIRanges <- function(gr,win,step,limit,maxWin){
  startWin <- ceiling((start(gr)-win+(1-limit)*width(gr))/step)+1
  startWin[startWin<1] <- 1
  startWin[startWin>maxWin] <- maxWin
  endWin <- floor((end(gr)-1-(1-limit)*width(gr))/step)+1
  endWin[endWin>maxWin] <- maxWin
  IRanges(start = startWin,end = endWin)
}
