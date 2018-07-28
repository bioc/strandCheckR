#' @title Get the Sliding Windows from a GRanges object
#' 
#' @description Get the positive/negative windows that overlap a GRanges object
#' 
#' 
#' @param x a GRanges object
#' @param seqInfo a data frame that contains some key information of 
#' the alignments
#' @param winWidth The width of each window
#' @param winStep The step size for sliding the window
#' @return A list of two logical vectors (for positive and negative strand) 
#' defining which windows that overlap the given Granges objects
#' 
#' @importFrom GenomicRanges start<-
#' @importFrom GenomicRanges end<-
#' @importFrom GenomicRanges ranges<-
#' @importFrom BiocGenerics strand
#' @importFrom GenomeInfoDb seqlevels
#' 

getWinFromGranges <- function(x, seqInfo, winWidth = 1000, winStep = 100) 
{   
    # Check the correct columns are in the seqInfo df
    reqCols <- c("FirstBaseInPart")
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

