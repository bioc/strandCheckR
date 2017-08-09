#' @title Partition the chromosomes
#' 
#' @description Partition the set of chromosomes for easier processing
#' 
#' @details Partition the set of chromosomes such that each part contains a 
#' subset of chromosome whose sum of length is >= \code{partitionSize}.
#' This is to enable faster processing by \code{getWinFromBamFile} and \code{filterDNA}
#' 
#' @param chromosomes a vector of chromosomes to be partitioned
#' @param lengthSeq the length of each chromosome
#' @param partitionSize the minimum sum of length of the chromosomes in each partition
#' @export
#'
partitionChromosomes <- function(chromosomes, lengthSeq, partitionSize){
  sumLength <- cumsum(as.numeric(lengthSeq))
  partition <- list()
  currentSum <- 0
  firstChrInPart <- 1
  lastChrInPart <- 1
  idPart <- 1
  while (lastChrInPart <= length(chromosomes)){
    lastChrInPart <- which(sumLength >= currentSum+partitionSize)
    if (length(lastChrInPart) == 0){
      lastChrInPart <- length(chromosomes)
    }
    else{
      lastChrInPart <- lastChrInPart[1]
    }
    partition[[idPart]] <- chromosomes[firstChrInPart:lastChrInPart]
    currentSum <- sumLength[lastChrInPart]
    firstChrInPart <- lastChrInPart+1
    lastChrInPart <- lastChrInPart + 1
    idPart <- idPart + 1
  }
  return(partition)
}
