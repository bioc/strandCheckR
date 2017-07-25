#' @title get window data frame with the correct chromosome name and position
#' 
#' @description get the correct chromsome name and position for each window
#' 
#' @param Win a data frame contains the information of every window in \code{chromosomes}
#' @param chromsomes a vector of chromosome names
#' @param statInfo a data frame that contains some information of the alignments
#' @param win the length of sliding window
#' @param step the step length to sliding the window
#' @export
#'
getWinInChromosome <- function(Win,chromosomes,statInfo,win,step){
  Chromosome <- rep("",nrow(Win))
  for (i in seq_along(chromosomes)){
    if (!is.na(statInfo$FirstBaseInPartition[i])){
      mi <- ceiling((statInfo$FirstBaseInPartition[i])/step)
      ma <- ceiling((statInfo$LastBaseInPartition[i]-win+1)/step)
      j <- which(Win$Start >=mi & Win$Start <=ma)
      Chromosome[j] <- chromosomes[i]
      Win$Start[j] <- Win$Start[j] - mi +1  
    }
  }
  Win[["Chr"]] <- Chromosome  
  return(Win)
}

