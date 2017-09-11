#' @title Plot the histogram of positive proportions of the input windows data frame
#'
#' @description Plot the histogram of positive proportions of the input windows data frame
#'
#' @param windows data frame containing the number of positive/negative reads for each window
#' Windows can be get by the function \code{getWinFromBamFile}.
#' @param breaks an integer vector that specifies how you want to partition the windows based on the number of reads. By default \code{breaks} = c(10,100,1000), which means that your windows will be partition into 4 groups, those have number of reads < 10, from 10 to 100, from 100 to 1000, and > 1000
#' @param save if TRUE, then the plot will be save into the file given by \code{file} parameter
#' @param file the file name to save to plot
#' @param facet_wrap_chromosomes if TRUE, then the plots will be splitted by chromosomes. FALSE by default
#' @param useCoverage if TRUE then plot the coverage strand information, otherwise plot the number of reads strand information. FALSE by default
#' @param ... used to pass parameters to facet_wrap
#' @seealso getWinFromBamFile, plotWin
#'
#' @examples
#' \dontrun{
#' #for single end bam file
#' bamfilein = system.file("extdata","s1.chr1.bam",package = "strandCheckR")
#' windows <- getWinFromBamFile(file = bamfilein)
#' plotHist(windows)
#' 
#' #for paired end bamfile
#' bamfilepair = system.file("extdata","120.bam",package = "strandCheckR")
#' windowsP <- getWinFromBamFile(file = bamfilepair)
#' plotHist(windowsP)
#' }
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom grid unit
#' @importFrom ggplot2 facet_wrap
#' @importFrom magrittr set_colnames
#' @importFrom dplyr mutate select filter one_of starts_with
#' @importFrom stringr str_extract
#' @export
plotHist <- function(windows, breaks=c(10,100,1000), save=FALSE, file = "hist.pdf", 
                     facet_wrap_chromosomes=FALSE, useCoverage=FALSE, ...){
  
  # The initial checks for appropriate input
  reqWinCols <- c("Chr", "Start", "NbPositive", "NbNegative", "CovPositive", "CovNegative", "MaxCoverage")
  stopifnot(all(reqWinCols %in% colnames(windows)))
  stopifnot(is.numeric(breaks))
  stopifnot(is.logical(c(save, facet_wrap_chromosomes, useCoverage)))
  
  # Check to see if we have SE or PE
  readType <- ifelse("Type" %in% colnames(windows), "PE", "SE")
  if (readType == "SE") windows$Type <- "R1"
  
  # Make sure we don't facet if we only have one chromosome
  nChr <- length(unique(windows$Chr))
  if (nChr == 1) facet_wrap_chromosomes <- FALSE
  
  # Calculate the proportion of reads for the + strand, based on either coverage or the number of reads
  keepCols <- c("Nb", "Cov")[useCoverage + 1]
  windows <- as.data.frame(windows)
  windows <- dplyr::select(windows, one_of(c("Chr", "MaxCoverage", "Type")), starts_with(keepCols))
  names(windows) <- str_extract(names(windows),"Chr|MaxCoverage|Type|Pos|Neg")
  windows$PositiveProportion <- windows$Pos / (windows$Pos + windows$Neg)

  
  # Annotate windows by the level of coverage in the each window
  if (length(breaks)==0){
    windows$group <- as.factor("all")
  } 
  else{
    covBreaks <- unique(c(0, breaks, max(windows$MaxCoverage)))
    nBreaks <- length(covBreaks)
    breakMat <- cbind(covBreaks[-nBreaks], covBreaks[-1])
    covLabels <- apply(breakMat, MARGIN = 1, FUN = function(x){paste0("(", x[1], ",", x[2], "]")})
    covLabels[nBreaks - 1] <- gsub("\\(([0-9]+),.+", "> \\1", covLabels[nBreaks-1])
    windows$group <- cut(windows$MaxCoverage, breaks = covBreaks, include.lowest = TRUE, labels = covLabels)
  }
  
  # Group the stranded proportions into 102 bins
  windows$propReads <- cut(windows$PositiveProportion, breaks = seq(0, 1, length.out = 102), include.lowest = TRUE)
  # Tally the bins by coverage
  windows <- reshape2::dcast(windows, Type + propReads ~ group, fun.aggregate = length, value.var = "propReads")
  # Melt for easy plotting with ggplot2
  windows <- reshape2::melt(windows, id.vars = c("Type", "propReads"), variable.name = "Coverage")
  windows$propReads <- (as.integer(windows$propReads) - 1) / 100
  # Convert the numbers of windows into proportions keeping R1 & R2 separate
  windows <- lapply(split(windows, f = windows$Type), function(x){
    x$value <- x$value / sum(x$value)
    x
  })
  windows <- dplyr::bind_rows(windows)
  
  # Organise the faceting
  facets <- c()
  if (readType == "PE" && facet_wrap_chromosomes){
    facets <- ~Chr + Type
  }
  if (readType == "PE" && !facet_wrap_chromosomes){
    facets <- ~Type
  }
  if (readType == "SE" && facet_wrap_chromosomes){
    facets <- ~Chr
  }
  
  # Make the plot
  g <- ggplot(windows, aes_string(x = "propReads", y = "value", fill = "Coverage")) + 
    geom_bar(stat = "identity") + 
    labs(x = "Proportion of Reads on '+' Strand",
         y = "Proportion of Windows") +
    theme_bw() +
    theme(plot.margin = unit(c(0.02,0.04,0.03,0.02), "npc"))
  if (!is.null(facets)) {
    
    # Get any facet arguments from dotArgs that have been set manually
    dotArgs <- list(...)
    allowed <- names(formals(ggplot2::facet_wrap))
    keepArgs <- names(dotArgs) %in% setdiff(allowed, "facets")
    argList <- c(list(facets = facets), dotArgs[keepArgs])
    myFacets <- do.call(facet_wrap, argList)
    g <- g + myFacets
    
  }
    
  if (save==TRUE){
    message("The plot will be saved to the file ",file)
    ggplot2::ggsave(filename = file,plot = g)
  }

  g

}

