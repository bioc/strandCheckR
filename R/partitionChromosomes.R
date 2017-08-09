#' @title partition the set of chromosomes such that each part contains a subset of chromosome whose sum of length is at least 1e8
#' @param chromosomes a vector of chromosomes to be partitionned
#' @param lengthSeq the length of each chromosome
#' @param yieldSize the minimum sum of length of the chromosomes in each partition
#' @export
#'
partitionChromosomes <- function(chromosomes,lengthSeq,yieldSize){
  sumLength <- cumsum(as.numeric(lengthSeq))
  partition <- list()
  currentSum <- 0
  firstChrInPart <- 1
  lastChrInPart <- 1
  idPart <- 1
  while (lastChrInPart <= length(chromosomes)){
    lastChrInPart <- which(sumLength >= currentSum+yieldSize)
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
