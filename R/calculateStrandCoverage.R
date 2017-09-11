#' @title Calculate the strand information based on coverage
#'
#' @description Calculate the coverage coming from '+'/'-' reads in all sliding wndows
#' 
#' @param winPositiveAlignments a list that has a `Coverage` field that contains coverage coming from positive reads
#' @param winNegativeAlignmentsa a list that has a `Coverage` field that contains coverage coming from negative reads
#' @param winWidth the length of the sliding window, 1000 by default.
#' @param winStep the step length to sliding the window, 100 by default.
#' @importFrom IRanges Views
#' @importFrom S4Vectors Rle


calculateStrandCoverage <- function(winPositiveAlignments, winNegativeAlignments, winWidth = 1000, winStep= 100){
  # make sure winPositiveAlignments$Coverage and winNegativeAlignments$Coverage
  # have the same length to avoid some warnings afterward
  lastBase <- max(c(length(winPositiveAlignments$Coverage),
                    length(winNegativeAlignments$Coverage)))
  lenPC <- length(winPositiveAlignments$Coverage)
  lenNC <- length(winNegativeAlignments$Coverage)
  if (lenNC < lastBase) {
    winNegativeAlignments$Coverage <- c(winNegativeAlignments$Coverage,rep(0, lastBase - lenNC))
  } 
  if (lenPC < lastBase){
    winPositiveAlignments$Coverage <- c(winPositiveAlignments$Coverage,rep(0, lastBase - lenPC))
  }
  
  #calculate the number of positive and negative bases in each window
  nbWin <- ceiling((lastBase - winWidth) / winStep) + 1
  st <- seq(1, (nbWin - 1)*winStep + 1, by =winStep)
  end <- seq(winWidth, (nbWin-1)*winStep + winWidth, by = winStep)
  CovPositive <- Views(winPositiveAlignments$Coverage, start = st, end = end) 
  CovPositive <- Rle(sum(CovPositive))
  CovNegative <- Views(winNegativeAlignments$Coverage, start = st, end = end) 
  CovNegative <- Rle(sum(CovNegative))
  
  #calculate max the max coverage in each window
  maxCoverage <- Rle(max(Views(winPositiveAlignments$Coverage + winNegativeAlignments$Coverage, 
                               start = st, end=end)))
  list("CovPositive"=CovPositive,"CovNegative"=CovNegative,"MaxCoverage"=maxCoverage)
  
}