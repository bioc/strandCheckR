#' @title get the kept probability of each read fragment
#' @description  calculate the keeping probability of each read fragment based on the keeping probability of the windows containing it
#' @param fragments an IRange object defind the starting, ending position of each fragment
#' @param keptProbaW an Rle object define the kept probability of every sliding window
#' @param errorRate the rate that a read has wrong strand
#' @export
#'

keptAlignment <- function(fragments,keptProbaW,errorRate){
    if(length(fragments)>0){
      minL <- min(width(fragments))
      maxL <- max(width(fragments))
      proba <- rep(0,length(fragments))
      for (m in minL:maxL){
        id <- which(width(fragments)==m)
        x <- keptProbaW[start(fragments)[id]]
        for (j in seq_along(1:m)[-m]){
          x <- myMax(x,keptProbaW[j+start(fragments)[id]])
        }
        suppressWarnings(proba[id] <- rbinom(length(id),1,as.vector(x)))
        rm(id)
        rm(x)
      }
      return(which(proba==1))
    }
    else{
      return(c())
    }
}

