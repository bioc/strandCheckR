#' @title Plot the number of reads vs positive proportion of the input windows data frame.
#'
#' @description Plot the number of reads vs positive proportion of the input windows data frame.

#' @param windows data frame containing the number of positive/negative reads for each window.
#' Windows can be get by the function \code{getWinFromBamFile}.
#'
#' @param breaks an integer vector that specifies how you want to partition the windows based on the number of reads. By default \code{breaks} = c(10,100,1000), which means that your windows will be parition into 4 groups, those have number of reads < 10, from 10 to 100, from 100 to 1000, and > 1000
#'
#' @param threshold a real vector that specifies which threshold lines you want to draw on the plot. The positive window above the threshold line (or negative window below the threshold line) will be kept if we use that threshold to filter the data. By default, the thresholds 0.6, 0.7, 0.8, 0.9 will be plotted.
#'
#' @param pvalue 0.05 by default, any window has pvalue below it in the test comparing with a threshold will be kept.
#' @param save if TRUE, then the plot will be save into the file given by \code{file} parameter
#' @param file the file name to save to plot
#' @param facet_wrap_chromosomes if TRUE, then the plots will be splitted by chromosomes. FALSE by default
#' @param useCoverage if TRUE then plot the coverage strand information, otherwise plot the number of reads strand information. FALSE by default
#' @seealso getWinFromBamFile,  plotHist
#'
#' @importFrom dplyr select mutate distinct one_of starts_with
#' @importFrom stats pnorm
#' 
#' @examples
#' \dontrun{
#' bamfilein = system.file("extdata","s1.chr1.bam",package = "strandCheckR")
#' windows <- getWinFromBamFile(file = bamfilein)
#' plotWin(windows)
#' @export
#'

plotWin <- function(windows,breaks=c(10,100,1000),threshold=c(0.6,0.7,0.8,0.9),pvalue=0.05,save=FALSE,file="win.pdf",facet_wrap_chromosomes=FALSE,useCoverage=FALSE){
  
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
  
  #Remove duplicated points for lightening the plots
  if (readType == "PE" && facet_wrap_chromosomes){
    windowsReduced <- windows %>% mutate(NbReads = round(NbReads, -1), PositiveProportion = round(PositiveProportion, 2)) %>%
      distinct(Chr,Type,NbReads,PositiveProportion,group) 
  }
  if (readType == "PE" && !facet_wrap_chromosomes){
    windowsReduced <- windows %>% mutate(NbReads = round(NbReads, -1), PositiveProportion = round(PositiveProportion, 2)) %>%
      distinct(Type,NbReads,PositiveProportion,group) 
  }
  if (readType == "SE" && facet_wrap_chromosomes){
    windowsReduced <- windows %>% mutate(NbReads = round(NbReads, -1), PositiveProportion = round(PositiveProportion, 2)) %>%
      distinct(Chr,NbReads,PositiveProportion,group) 
  }
  if (readType == "SE" && !facet_wrap_chromosomes){
    windowsReduced <- windows %>% mutate(NbReads = round(NbReads, -1), PositiveProportion = round(PositiveProportion, 2)) %>%
      distinct(NbReads,PositiveProportion,group) 
  }
  rm(windows)
  
  # Calculating the threshold lines
  ThresholdP <- data.frame("NbReads" = c(), "PositiveProportion" = c(), "Threshold"= c())
  ThresholdN <- data.frame("NbReads" = c(), "PositiveProportion" = c(), "Threshold"= c())
  nbSampling <- 10000
  for (t in threshold){
    positiveReadsT <- sapply(1:nbSampling,function(N){
      tP = log(t/(1-t))
      p = seq(round(N*t),N,1)
      pP = p/N
      mP = log(pP/(1-pP))
      sdP = sqrt(1/(N*pP*(1-pP)))
      pNorm <- pnorm(tP,mean = mP, sd = sdP)#pBinom <- pbinom(t*N,size = N, prob = pP)
      aNorm <- which(pNorm <= 0.05)[1]
      return(p[aNorm])
    })
    tP <- data.frame("NbReads" = 1:nbSampling, "PositiveProportion" = positiveReadsT/(1:nbSampling), "Threshold"= paste0(t)) #%>% rbind(data.frame("NbPositive"=max(windowsReduced$NbPositive),"NbNegative"=round(max(windowsReduced $NbPositive)/t)-max(windowsReduced $NbPositive),"Threshold"=paste0(t)))
    ThresholdP <- rbind(ThresholdP,tP)
    tN <- data.frame("NbReads" = 1:nbSampling, "PositiveProportion" = 1-positiveReadsT/(1:nbSampling), "Threshold"= paste0(t)) #%>% rbind(data.frame("NbPositive"=round(max(windowsReduced $NbNegative)/t)-max(windowsReduced $NbNegative),"NbNegative"=max(windowsReduced $NbNegative),"Threshold"=paste0(t)))
    ThresholdN <- rbind(ThresholdN,tN)
  }
  if (max(windowsReduced$NbReads) > nbSampling){
    ThresholdP <- rbind(ThresholdP,data.frame("NbReads" = max(windowsReduced$NbReads), "PositiveProportion" = threshold, "Threshold"=paste0(threshold)))
    ThresholdN <- rbind(ThresholdN,data.frame("NbReads" = max(windowsReduced$NbReads), "PositiveProportion" = (1-threshold), "Threshold"=paste0(threshold)))
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
  
  g <- ggplot2::ggplot() +
    ggplot2::geom_point(data = windowsReduced, ggplot2::aes(x = NbReads, y = PositiveProportion, colour = group)) +
    ggplot2::geom_line(data = ThresholdP, ggplot2::aes(x = NbReads, y = PositiveProportion, linetype = Threshold)) +
    ggplot2::geom_line(data = ThresholdN, ggplot2::aes(x = NbReads, y = PositiveProportion, linetype = Threshold)) +
    labs(x = "Number of Reads",
         y = "Proportion of Reads on '+' Strand",
         colour = "Max Coverage") +
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
