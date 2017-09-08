#' @title Get the Sliding Windows from a GRanges object
#' 
#' @description Get the positive/negative windows that overlap a GRanges object
#' 
#' 
#' @param x a GRanges object
#' @param chromosomeInfo a data frame that contains some key information of the alignments
#' @param winWidth The width of each window
#' @param winStep The step size for sliding the window
#' 
#' @return A list of two logical vectors (for positive and negative strand) defining which windows that overlap the given Granges objects
#' 
#' @importFrom GenomicRanges start<-
#' @importFrom GenomicRanges end<-
#' @importFrom GenomicRanges ranges<-
#' @importFrom BiocGenerics strand
#' @importFrom GenomeInfoDb seqlevels
#' 
#' @export
getWinFromGranges <- function(x, chromosomeInfo, winWidth = 1000, winStep = 100){
  
  # Check for correct chromosomes 
  if (!all(chromosomes %in% seqlevels(x))) stop("Invalid specification of chromosomes.\nMust be in the set", seqlevels(x))
  # Handle missing chromsosomeInfo
  if (missing(chromosomes)) chromosomes <- seqlevels(x)
  # Check the correct columns are in the chromosomeInfo df
  reqCols <- c("FirstBaseInPartition")
  if (!all(reqCols %in% names(chromosomeInfo))) stop("chromosomeInfo must contain the column ", reqCols)
  stopifnot(is.numeric(winWidth) || is.numeric(winStep))
  
  for (i in seq_along(chromosomeInfo$Sequence)){
    r <- which(as.vector(seqnames(x)) == chromosomeInfo$Sequence[i])
    if (length(r)>0){
      start(ranges(x)[r]) <- start(ranges(x)[r]) + chromosomeInfo$FirstBaseInPartition[i] -1
      end(ranges(x)[r]) <- end(ranges(x)[r]) + chromosomeInfo$FirstBaseInPartition[i] -1  
    }
  }
  
  mustKeepPos <- getWinFromIRanges(ranges(x)[strand(x)!="-",], winWidth, winStep, 1) 
  mustKeepPos <- coverage(mustKeepPos) > 0
  mustKeepNeg <- getWinFromIRanges(ranges(x)[strand(x)!="+",], winWidth, winStep, 1) 
  mustKeepNeg <- coverage(mustKeepNeg) > 0
  
  list(Positive = mustKeepPos, Negative = mustKeepNeg)
}

