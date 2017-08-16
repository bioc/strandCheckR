calculateStrandNbReads <- function(winPositiveAlignments, winNegativeAlignments){
  
  # Calculate strand information based on number of reads 
  # have the same length to avoid some warnings afterward
  NbPositive <- coverage(winPositiveAlignments$Win)
  NbNegative <- coverage(winNegativeAlignments$Win)
  
  # Find the last window in both sets of windows
  lastWin <- max(c(end(winNegativeAlignments$Win), 
                   end(winPositiveAlignments$Win)))
  
  # Fill with zeroes if required
  #make sure NbPositive and NbNegative have the same length
  lenP <- length(NbPositive)
  lenN <- length(NbNegative)
  if (lenN < lastWin) NbNegative <- c(NbNegative,rep(0,lastWin - lenN))
  if (lenP < lastWin) NbPositive <- c(NbPositive,rep(0,lastWin - lenP)) 
  
  list("NbPositive"=NbPositive,"NbNegative"=NbNegative)
}