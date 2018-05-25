#' @title Plot the histogram of positive proportions of the input histogram data frame
#'
#' @description Plot the histogram of positive proportions of the input histogram  data frame
#'
#' @param histWin data frame obtained from \code{summarizeHist} function
#' @param save if TRUE, then the plot will be save into the file given by \code{file} parameter
#' @param file the file name to save to plot
#' @param facets colnames of \code{histWin} which will be used to split the plot
#' @param heatmap if TRUE, then use heat map to plot the histogram, otherwise use barplot. FALSE by default.
#' @param ... used to pass parameters to facet_wrap
#' @seealso summarizeHist, getWinFromBamFile, plotWin
#'
#' @examples
#' \dontrun{
#' #for single end bam file
#' bamfilein = system.file("extdata","s1.sorted.bam",package = "strandCheckR")
#' windows  <- getWinFromBamFile(file = bamfilein)
#' histWin <- summarizeHist(windows)
#' plotHist(histWin)}
#' 
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_bar geom_tile
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom grid unit
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 ggsave ggtitle scale_fill_gradient
#' @importFrom magrittr set_colnames
#' @importFrom dplyr mutate select filter one_of starts_with
#' @importFrom stringr str_extract
#' @export
plotHist <- function(histWin, save=FALSE, file = "hist.pdf", facets = NULL, heatmap = FALSE, ...){
  
  # The initial checks for appropriate input
  reqWinCols <- c("PosStrandProp","ReadCountProp")
  stopifnot(all(reqWinCols %in% colnames(histWin)))
  stopifnot(is.logical(save))
  allows_facet_wrap <- setdiff(colnames(histWin),c(reqWinCols,"Coverage"))
  facets <- intersect(facets,allows_facet_wrap) 
  
  # Make the plot
  if (heatmap==FALSE){
    g <- ggplot(histWin, aes_string(x = "PosStrandProp", y = "ReadCountProp", fill = "Coverage")) + 
      geom_bar(stat = "identity") + 
      labs(x = "Proportion of Reads on '+' Strand",
           y = "Proportion of Windows") +
      theme_bw() +
      theme(plot.margin = unit(c(0.02,0.04,0.03,0.02), "npc"))  
    if (length(facets)>0) {
      
      # Get any facet arguments from dotArgs that have been set manually
      dotArgs <- list(...)
      allowed <- names(formals(facet_wrap))
      keepArgs <- names(dotArgs) %in% setdiff(allowed, "facets")
      argList <- c(list(facets = facets), dotArgs[keepArgs])
      myFacets <- do.call(facet_wrap, argList)
      g <- g + myFacets
    }
  } else{
    cov <- unique(histWin$Coverage)
    p <- list()
    if (!("File" %in% colnames(histWin))){
      histWin$File <- "File"
    }
    facets <- facets[facets!="File"]
    if (length(facets)>0) {
      # Get any facet arguments from dotArgs that have been set manually
      dotArgs <- list(...)
      allowed <- names(formals(facet_wrap))
      keepArgs <- names(dotArgs) %in% setdiff(allowed, "facets")
      argList <- c(list(facets = facets), dotArgs[keepArgs])
      myFacets <- do.call(facet_wrap, argList)
    }
    for (i in 1:length(cov)){
      l <- filter(histWin,Coverage==cov[i]) 
      p[[i]] <- ggplot(l,aes_string(x="PosStrandProp",y="File",fill="ReadCountProp")) + 
        geom_tile() +
        scale_fill_gradient(low="white",high="red",na.value = "white")
      if (length(facets)>0){
        p[[i]] <- p[[i]] + myFacets
      }
      if (length(cov)>1) {
        p[[i]] <- p[[i]] + ggtitle(paste0("Coverage ",cov[i]))
      }
    }
    g <- grid.arrange(grobs=p, nrow = ceiling(length(cov)/2))
  }
  
  if (save==TRUE){
    message("The plot will be saved to the file ",file)
    dotArgs <- list(...)
    allowed <- names(formals(ggsave))
    keepArgs <- names(dotArgs) %in% setdiff(allowed, c("filename","plot"))
    argList <- c(list(filename = file, plot = g), dotArgs[keepArgs])
    do.call(ggsave,argList)
  }
  if (heatmap==FALSE) g
}

