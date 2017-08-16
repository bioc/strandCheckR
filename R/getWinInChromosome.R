# #' @title get window data frame with the correct chromosome name and position
# #' 
# #' @description get the correct chromsome name and position for each window
# #' 
# #' @param Win a data frame contains the information of every window in \code{chromosomes}
# #' @param chromosomes a vector of chromosome names
# #' @param chromosomeInfo a data frame that contains some information of the alignments
# #' @param winWidth the length of sliding window
# #' @param winStep the step length to sliding the window
# #' 
getWinInChromosome <- function(Win, chromosomes, chromosomeInfo, winWidth = 1000, winStep = 100){
  
  # Check the correct columns are in the chromosomeInfo df
  reqCols <- c("FirstBaseInPartition", "LastBaseInPartition")
  if (!all(reqCols %in% names(chromosomeInfo))) stop("chromosomeInfo must contain the columns ", reqCols)
  stopifnot(is.numeric(winWidth) || is.numeric(winStep))
  
  for (i in seq_along(chromosomes)){
    if (!is.na(chromosomeInfo$FirstBaseInPartition[i])){
      currentChr <- chromosomes[i]
      id <- which(chromosomes == currentChr)
      idFirst <- ceiling(chromosomeInfo$FirstBaseInPartition[id] / winStep) # id of the first window of the chromosome
      idLast <- ceiling((chromosomeInfo$LastBaseInPartition[id] - winWidth+1) / winStep)# id of the last window of the chromosome
      idRows <- which(Win$Start >= idFirst & Win$Start <= idLast) #get the windows of the chromosome
      Win$Chr[idRows] <- Rle(currentChr)
      Win$Start[idRows] <- (Win$Start[idRows] - idFirst)*winStep +1
    }
  }
  
  return(Win)
}
