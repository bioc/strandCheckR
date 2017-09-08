#' @title Get the Ranges of Sliding Windows from an IRanges object
#' 
#' @description Get the Ranges of Sliding Windows from an IRanges object
#' 
#' @details
#' This finds the windows that overlap each fragment of a read and returns a range containing
#' this list of windows for each read fragment.
#' This allows the total number of read fragments within a window to be calculated simply using \link{coverage}.
#' 
#' @param x an IRanges object containing the start and end position of each fragment
#' @param winWidth The width of each window
#' @param winStep The step size for sliding the window
#' @param readProp A read is considered to be included in a window if at least \code{readProp} of it is in the window. 
#' Specified as a proportion.
#' @param maxWin The maximum window ID
#' 
#' @return An IRanges object containing the index of the windows containing each fragment
#' 
#' @export
#'
#'
getWinFromIRanges <- function(x, winWidth = 1000L, winStep = 100L, readProp = 0, maxWin = Inf){
  
  readProp <- readProp[1]
  if (readProp < 0 || readProp > 1) stop("readProp must be within the range [0, 1]")
  
  stopifnot(is.numeric(maxWin))
  
  # Calculate the index of the first window that overlaps a fragment
  startWin <- ceiling((start(x) - winWidth + (readProp)*width(x)) / winStep) + 1
  startWin[startWin < 1] <- 1
  startWin[startWin > maxWin] <- maxWin
  
  # Calculate the index of the final window that overlaps a fragment
  endWin <- floor((end(x) - 1 - (readProp)*width(x)) / winStep) + 1
  endWin[endWin > maxWin] <- maxWin
  # Return the complete set of windows for each fragment
  IRanges(start = startWin, end = endWin)
}
