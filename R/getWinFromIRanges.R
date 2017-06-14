#' @title get the ranges of windows from on IRanges
#'
#' @export
#'
#'
getWinFromIRanges <- function(gr,win,step,limit){
  startWin <- ceiling((start(gr)-win+(1-limit)*width(gr))/step)+1
  startWin[startWin<1] <- 1
  endWin <- floor((end(gr)-1-(1-limit)*width(gr))/step)+1
  IRanges(start = startWin,end = endWin)
}
