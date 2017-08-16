#' @title Get the Ranges of Sliding Windows from a GRanges object
#' 
#' @description Get the Ranges of positive/negative windows that overlap a GRanges object
#' 
#' 
#' @param x a GRanges object
#' @param chromosomes a list of chromosome names
#' @param chromosomeInfo a data frame that contains some information of the alignments
#' @param winWidth The width of each window
#' @param winStep The step size for sliding the window
#' 
#' @return A list of two logical vectors defining whether to keep windows based on positive or negative strands
#' 
#' @importFrom GenomicRanges start<-
#' @importFrom GenomicRanges end<-
#' @importFrom GenomicRanges ranges<-
#' @importFrom BiocGenerics strand
#' @importFrom GenomeInfoDb seqlevels
#' 
#' @export
getWinFromGranges <- function(x, chromosomes, chromosomeInfo, winWidth = 1000, winStep = 100){
  
  # Check for correct chromosomes 
  if (!all(chromosomes %in% seqlevels(x))) stop("Invalid specification of chromosomes.\nMust be in the set", seqlevels(x))
  # Handle missing chromsosomeInfo
  if (missing(chromosomes)) chromosomes <- seqlevels(x)
  # Check the correct columns are in the chromosomeInfo df
  reqCols <- c("FirstBaseInPartition")
  if (!all(reqCols %in% names(chromosomeInfo))) stop("chromosomeInfo must contain the column ", reqCols)
  stopifnot(is.numeric(winWidth) || is.numeric(winStep))
  
  for (i in seq_along(chromosomes)){
    r <- which(as.vector(seqnames(x)) == chromosomes[i])
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

