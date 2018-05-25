#' @title Intersect the windows with an annotation data frame to get feature for each window
#'
#' @description Intersect the windows with an annotation data frame to get feature for each window
#'
#' @param windows data frame containing the strand information of the sliding windows.
#' Windows can be obtained using the function \code{getWinFromBamFile}.
#' @param annotation a Grange object that you want to intersect with your windows. 
#' It can have mcols which contains the information or features that could be able to integrate into the input windows 
#' @param getFeatureInfo whether to get the information of features in the mcols of annotation data. 
#' If FALSE the return windows will have an additional column OverlapAnnot indicating wheather each window overlaps with any range in the annotion data.
#' If TRUE the return windows will contain the information of features that overlap each window
#' @param mcolsAnnot the column names of the mcols of annotation that you want to get information
#' @param collapse character which is used collapse multiple features that overlap with a same window into a string.
#' If missing then we don't collapse them.
#' @param ... used to pass parameters to GenomicRanges::findOverlaps
#' @seealso getWinFromBamFile, summarizeHist, plotHist, plotWin
#'
#' @examples
#' \dontrun{
#' #for single end bam file
#' bamfilein = system.file("extdata","s1.sorted.bam",package = "strandCheckR")
#' windows <- getWinFromBamFile(file = bamfilein)
#' windows$Chr <- "chr10"
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' annot <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' windows <- intersectWithFeature(windows,annot)
#' h <- summarizeHist(windows,group_by = "OverlapFeature")
#' plotHist(h,facets = "OverlapFeature")
#' plotWin(h,facets = "OveralapFeature")
#' }
#' 
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges GRanges mcols findOverlaps
#' @export


intersectWithFeature <- function(windows, annotation, getFeatureInfo = FALSE, overlapCol = "OverlapFeature", mcolsAnnot, collapse, ...){
  #check annotation is a GRanges object
  stopifnot(class(annotation)=="GRanges")
  #check if windows contains required column
  reqWinCols <- c("Chr","Start", "End")
  stopifnot(all(reqWinCols %in% colnames(windows)))
  w <- GRanges(seqnames = windows$Chr, ranges = IRanges(start = windows$Start, end = windows$End))
  if (length(intersect(seqlevels(w),seqlevels(annotation))) == 0){
    message("Windows do not have any sequence in annotation data.")
    return(windows)
  }
  
  #Get overlap between windows and annotation
  ol <- findOverlaps(w,annotation, select = "all", ...)
  
  windows$OverlapFeature <- "NotOverlap"
  windows$OverlapFeature[as.integer(unique(from(ol)))] <- "Overlap"
  if (getFeatureInfo){
    
    ## check collapse
    stopifnot(!missing(collapse) && is.character(collapse))
    
    ## Check mcolsAnnot
    if (missing(mcolsAnnot)){
      mcolsAnnot <- mcols(annotation)
    } else{
      if (length(mcolsAnnot)==0){
        message("No column specified to get from mcols of annotation.")
        return(windows)
      } else{
        mcolsAnnot <- intersect(mcolsAnnot,colnames(mcols(annotation)))    
        if (length(mcolsAnnot)==0){
          message("mcols of annotation does not contain any column ",mcolsAnnot)
          return(windows)
        }
      }
    }
    
    featureTo <- mcols(annotation)[to(ol),mcolsAnnot[1]]
    splitTo <- split(featureTo, f = from(ol), drop = TRUE)
    if (!missing(collapse)){
      splitTo <- sapply(splitTo,function(s){
        paste0(unique(s),collapse = collapse)
      })
    }
    indexFrom <- as.integer(names(splitTo))
    windows[[mcolsAnnot[1]]] <- "unknown"
    windows[[mcolsAnnot[1]]][indexFrom] <- splitTo 
    
    for (i in seq_along(mcolsAnnot)[-1]){
      featureTo <- mcols(annotation)[to(ol),mcolsAnnot[i]]
      splitTo <- split(featureTo, f = from(ol))
      if (!missing(collapse)){
        splitTo <- sapply(splitTo,function(s){
          paste0(unique(s),collapse = collapse)
        })
      }
      windows[[mcolsAnnot[i]]] <- "unknown"
      windows[[mcolsAnnot[i]]][indexFrom] <- splitTo 
    }
  }
  return(windows)
}