#' @title Get the Sliding Windows that overlap a GRanges object
#' 
#' @description Get the sliding windows that overlap a GRanges object. 
#' 
#' @details 
#' This finds the windows that overlaps the positive/negative strand of a 
#' GRanges object. The GRanges object, which is \code{mustKeepRanges} in the 
#' \code{filterDNA} method, defines the coordinates of the ranges in the 
#' reference genome that all reads mapped to those ranges must be kept in the 
#' filtering method. 
#' This method makes use of the \code{getWinFromIRanges} method by pretending 
#' the given range as the range of a read, and giving readProp=0 because we 
#' want to keep the reads as much as possible - so every read that overlaps 
#' even by one base with the given \code{mustKeepRanges}.
#' 
#' @param x a GRanges object, which defines the coordinates of 
#' the ranges in the reference genome that all reads mapped to those ranges
#' must be kept in the filtering method.
#' @param seqInfo a data frame that contains some key information of 
#' the alignments
#' @param winWidth the width of the sliding window, 1000 by default.
#' @param winStep the step length to sliding the window, 100 by default.
#' @return A list of two logical vectors (for positive and negative strand) 
#' defining which windows that overlap the given Granges objects
#' 
#' @importFrom GenomicRanges start<-
#' @importFrom GenomicRanges end<-
#' @importFrom GenomicRanges ranges<-
#' @importFrom BiocGenerics strand
#' @importFrom GenomeInfoDb seqlevels
#' 

getWinFromGRanges <- function(
    x, seqInfo, winWidth = 1000L, winStep = 100L
) 
{   
    # Check the correct columns are in the seqInfo df
    reqCols <- c("Sequence","FirstBaseInPart")
    if (!all(reqCols %in% names(seqInfo))) {
        stop("seqInfo must contain the column ", reqCols)
    }
    stopifnot(is.numeric(winWidth) || is.numeric(winStep))
    
    # Calculate start/end position of each sequence in the partition
    for (i in seq_along(seqInfo$Sequence)) {
        r <- which(as.vector(seqnames(x)) == seqInfo$Sequence[i])
        if (length(r) > 0) {
            start(ranges(x)[r]) <- 
                start(ranges(x)[r]) + seqInfo$FirstBaseInPart[i] - 1
            end(ranges(x)[r]) <- 
                end(ranges(x)[r]) + seqInfo$FirstBaseInPart[i] - 1
        }
    }
    
    # Calculate the windows that overlap the '+' ranges of x
    mustKeepPos <- getWinFromIRanges(
        ranges(x)[strand(x) != "-", ], winWidth, winStep, 0
        )
    mustKeepPos <- coverage(mustKeepPos) > 0
    
    # Calculate the windows that overlap the '-' ranges of x
    mustKeepNeg <- getWinFromIRanges(
        ranges(x)[strand(x) != "+", ], winWidth, winStep, 0
        )
    mustKeepNeg <- coverage(mustKeepNeg) > 0
    
    list(Positive = mustKeepPos, Negative = mustKeepNeg)
}

