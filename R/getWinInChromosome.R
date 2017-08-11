#' @title get window data frame with the correct chromosome name and position
#' 
#' @description get the correct chromsome name and position for each window
#' 
#' @param Win a data frame contains the information of every window in \code{chromosomes}
#' @param chromosomes a vector of chromosome names
#' @param statInfo a data frame that contains some information of the alignments
#' @param winWidth the length of sliding window
#' @param winStep the step length to sliding the window
#' @export
#'
getWinInChromosome <- function(Win,chromosomes,statInfo,winWidth,winStep){
  Win$Chr <- rep("",nrow(Win))
  for (i in seq_along(chromosomes)){
    if (!is.na(statInfo$FirstBaseInPartition[i])){
      currentChr <- chromosomes[i]
      id <- which(chromosomes == currentChr)
      idFirst <- ceiling(statInfo$FirstBaseInPartition[id] / winStep) # id of the first window of the chromosome
      idLast <- ceiling((statInfo$LastBaseInPartition[id] - winWidth+1) / winStep)# id of the last window of the chromosome
      idRows <- which(Win$Start >= idFirst & Win$Start <= idLast) #get the windows of the chromosome
      Win$Chr[idRows] <- currentChr
      Win$Start[idRows] <- (Win$Start[idRows] - idFirst)*winStep +1
    }
  }
  return(Win)
}
