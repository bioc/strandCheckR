calculateStrandCoverage <- function(winPositiveAlignments,winNegativeAlignments,winWidth,winStep){
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