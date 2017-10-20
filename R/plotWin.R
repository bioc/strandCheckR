#' @title Plot the Proportion of + Stranded Reads
#'
#' @description Plot the number of reads vs the proportion of '+' stranded reads.
#' 
#' @details This function will plot the proportion of '+' stranded reads for each window,
#'  against the number of reads in each window.
#'  The threshold lines indicate the hypothetical boundary where windows will contain reads to kept or discarded
#'  using the filtering methods of \code{\link{filterDNA}}.
#'  
#'  Any plot can be easily modified using standard ggplot2 syntax (see Examples).
#'  
#'  @return The plot will be returned as a standard ggplot2 object
#' 
#' @param windows data frame containing the strand information of the sliding windows.
#' Windows should be obtained using the function \code{\link{getWinFromBamFile}} to ensure the correct data structure.
#' @param breaks an integer vector that specifies how you want to partition the windows based on coverage. 
#' By default \code{breaks} = c(10,100,1000), partition windows into 4 groups based on these values.
#' @param threshold a \code{numeric} vector between 0.5 & 1 that specifies which threshold lines to draw on the plot. 
#' The positive windows above the threshold line (or negative windows below the threshold line) will be kept when using \code{\link{filterDNA}}. 
#' @param save if TRUE, then the plot will be save into the file given by \code{file} parameter
#' @param file the file name to save to plot
#' @param facet_wrap_chromosomes if TRUE, then the plots will be splitted by chromosomes. FALSE by default
#' @param useCoverage if TRUE then plot the coverage strand information, otherwise plot the number of reads strand information. FALSE by default
#' @param ... used to pass parameters to facet_wrap during plotting
#' @seealso \code{\link{getWinFromBamFile}},  \code{\link{plotHist}}
#'
#' @importFrom dplyr select mutate distinct one_of starts_with
#' @importFrom stats pnorm
#' @importFrom stringr str_extract
#' @importFrom ggplot2 ggplot geom_point aes_string labs theme_bw theme facet_wrap geom_line ggsave
#' @importFrom grid unit
#' @examples
#' \dontrun{
#' bamfilein = system.file("extdata","s1.chr1.bam",package = "strandCheckR")
#' windows <- getWinFromBamFile(file = bamfilein)
#' plotWin(windows)
#' 
#' # Change point colour using ggplot2
#' library(ggplot2)
#' plotWin(bamWindows) + scale_colour_manual(values = rgb(seq(0, 1, length.out = 4), 0, 0))}
#' 
#' @export
#'
plotWin <- function(windows, breaks=c(10,100,1000), threshold=c(0.6,0.7,0.8,0.9), 
                    save=FALSE, file="win.pdf", facet_wrap_chromosomes=FALSE, useCoverage=FALSE, ...){
  
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
  windows$NbReads <- windows$Pos + windows$Neg
  windows$PositiveProportion <- windows$Pos / windows$NbReads
  
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
  
  # Remove duplicated points for faster plotting
  windows <- mutate(windows, NbReads = round(NbReads, -1), PositiveProportion = round(PositiveProportion, 2))
  if (facet_wrap_chromosomes){
    windows <- distinct(windows, Chr, Type, NbReads, PositiveProportion, group) 
  }
  else {
    windows <- distinct(windows, Type, NbReads, PositiveProportion, group) 
  }
  
  # Generating the threshold lines
  maxReads <- max(windows$NbReads)
  ThresholdP <- data.frame("NbReads" = c(), "PositiveProportion" = c(), "Threshold"= c())
  ThresholdN <- data.frame("NbReads" = c(), "PositiveProportion" = c(), "Threshold"= c())
  nbSampling <- 1000
  x <- floor(seq(1, sqrt(maxReads), length.out = nbSampling)^2)
  for (t in threshold){
    tP = log(t/(1-t))
    positiveReadsT <- sapply(x,function(N){
      p = seq(round(N*t),N,1) # Number of positive reads
      pP = p/N # Proportion of positive reads
      mP = log(pP/(1-pP)) # Mean prop-pos-reads (logit scale)
      sdP = sqrt(1/(N*pP*(1-pP))) # SD prop-pos-reads
      pNorm <- pnorm(tP,mean = mP, sd = sdP)#pBinom <- pbinom(t*N,size = N, prob = pP)
      aNorm <- which(pNorm <= 0.05)[1]
      return(p[aNorm])
    })
    tP <- data.frame("NbReads" = x, "PositiveProportion" = positiveReadsT/(x), "Threshold"= paste0(t)) #%>% rbind(data.frame("NbPositive"=max(windows$NbPositive),"NbNegative"=round(max(windows $NbPositive)/t)-max(windows $NbPositive),"Threshold"=paste0(t)))
    ThresholdP <- rbind(ThresholdP,tP)
    tN <- data.frame("NbReads" = x, "PositiveProportion" = 1-positiveReadsT/(x), "Threshold"= paste0(t)) #%>% rbind(data.frame("NbPositive"=round(max(windows $NbNegative)/t)-max(windows $NbNegative),"NbNegative"=max(windows $NbNegative),"Threshold"=paste0(t)))
    ThresholdN <- rbind(ThresholdN,tN)
  }
  
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
  
  # Make the basic plot
  xlab <- c("Number of reads", "Number of aligned bases")[useCoverage + 1]
  g <- ggplot() +
    geom_point(data = windows, aes_string(x = "NbReads", y = "PositiveProportion", colour = "group")) +
    geom_line(data = ThresholdP, aes_string(x = "NbReads", y = "PositiveProportion", linetype = "Threshold")) +
    geom_line(data = ThresholdN, aes_string(x = "NbReads", y = "PositiveProportion", linetype = "Threshold")) +
    labs(x = xlab,
         y = "Proportion of Reads on '+' Strand", 
         colour = "Max Coverage") + 
    theme_bw() +
    theme(plot.margin = unit(c(0.02,0.04,0.03,0.02), "npc"))
  
  if (!is.null(facets)) {
    # Get any facet arguments from dotArgs that have been set manually
    dotArgs <- list(...)
    allowed <- names(formals(facet_wrap))
    keepArgs <- names(dotArgs) %in% setdiff(allowed, "facets")
    argList <- c(list(facets = facets), dotArgs[keepArgs])
    myFacets <- do.call(facet_wrap, argList)
    g <- g + myFacets
  }
  if (save==TRUE){
    message("The plot will be saved to the file ",file)
    ggsave(filename = file,plot = g)
  }
  g
}
