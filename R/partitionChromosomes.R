#' @title partition the set of chromosomes such that each part contains a subset of chromosome whose sum of length is at least 1e8
#'
#' @export
#'
partitionChromosomes <- function(chromosomes,lengthSeq,yieldSize){
  sumLength <- cumsum(as.numeric(lengthSeq))
  partition <- list()
  currentSum <- 0
  old_id <- 1
  i <- 1
  id <- 1
  while (id<=length(lengthSeq)){
    id <- which(sumLength>=currentSum+yieldSize)
    if (length(id)==0){
      id <- length(lengthSeq)
    }
    else{
      id <- id[1]
    }
    partition[[i]] <- chromosomes[old_id:id]
    old_id <- id+1
    currentSum <- sumLength[id]
    i <- i + 1
    id <- id + 1
  }
  return(partition)
}
