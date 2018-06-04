#' @title Partition the sequences
#' 
#' @description Partition the set of sequences for easier processing
#' 
#' @details Partition the set of sequences such that each part contains a 
#' subset of sequences whose sum of length is >= \code{partitionSize}.
#' This is to enable faster processing by \code{getWinFromBamFile} and 
#' \code{filterDNA}
#' 
#' @param sequences a vector of sequences to be partitioned
#' @param sq a seqinfo object to be partitioned
#' @param lengthSeq the length of each sequence
#' @param partitionSize the minimum sum of length of the sequences in 
#' each partition
#' 
#' @return a list whose each element is a character vector containing the 
#' sequence names of a part
#' 

partitionSequences <- function(sequences, lengthSeq, partitionSize = 1e8){

    stopifnot(length(sequences) == length(lengthSeq))
    sumLength <- cumsum(as.numeric(lengthSeq))
    partition <- list()
    currentSum <- 0
    firstSeqInPart <- 1
    lastSeqInPart <- 1
    idPart <- 1
    # Break into partitions such that the minimum size is > partitionSize
    while (lastSeqInPart <= length(sequences)){
        lastSeqInPart <- which(sumLength >= currentSum+partitionSize)
        if (length(lastSeqInPart) == 0){
            lastSeqInPart <- length(sequences)
        }
        else{
            lastSeqInPart <- lastSeqInPart[1]
        }
        partition[[idPart]] <- sequences[firstSeqInPart:lastSeqInPart]
        currentSum <- sumLength[lastSeqInPart]
        firstSeqInPart <- lastSeqInPart+1
        lastSeqInPart <- lastSeqInPart + 1
        idPart <- idPart + 1
    }
    return(partition)
}



#' @rdname partitionSequences
partitionSeqinfo <- function(sq, partitionSize = 1e8){

    stopifnot(class(sq) == "Seqinfo")
    lengthSeq <- seqlengths(sq)
    sequences <- seqnames(sq)
    sumLength <- cumsum(as.numeric(lengthSeq))
    partition <- list()
    currentSum <- 0
    firstSeqInPart <- 1
    lastSeqInPart <- 1
    idPart <- 1
    # Break into partitions such that the minimum size is > partitionSize
    while (lastSeqInPart <= length(sequences)){
        lastSeqInPart <- which(sumLength >= currentSum+partitionSize)
        if (length(lastSeqInPart) == 0){
            lastSeqInPart <- length(sequences)
        } else{
            lastSeqInPart <- lastSeqInPart[1]
        }
        partition[[idPart]] <- sequences[firstSeqInPart:lastSeqInPart]
        currentSum <- sumLength[lastSeqInPart]
        firstSeqInPart <- lastSeqInPart+1
        lastSeqInPart <- lastSeqInPart + 1
        idPart <- idPart + 1
    }
    return(partition)
}