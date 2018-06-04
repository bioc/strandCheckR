#' @title Plot the Proportion of + Stranded Reads
#'
#' @description Plot the number of reads vs the proportion of '+'
#'  stranded reads.
#' 
#' @details This function will plot the proportion of '+' stranded reads for 
#' each window,
#'  against the number of reads in each window.
#'  The threshold lines indicate the hypothetical boundary where windows will 
#'  contain reads to kept or discarded
#'  using the filtering methods of \code{\link{filterDNA}}.
#'  
#'  Any plot can be easily modified using standard ggplot2 syntax (see Examples)
#'  
#' @param windows data frame containing the strand information of the sliding 
#' windows. Windows should be obtained using the function 
#' \code{\link{getWinFromBamFile}} to ensure the correct data structure.
#' @param split an integer vector that specifies how you want to partition the 
#' windows based on coverage. 
#' By default \code{split} = c(10,100,1000), partition windows into 4 groups 
#' based on these values.
#' @param threshold a \code{numeric} vector between 0.5 & 1 that specifies 
#' which threshold lines to draw on the plot. 
#' The positive windows above the threshold line (or negative windows below the 
#' threshold line) will be kept when using \code{\link{filterDNA}}. 
#' @param save if TRUE, then the plot will be save into the file given by 
#' \code{file} parameter
#' @param file the file name to save to plot
#' @param facets colnames of \code{windows} which will be used to split the plot
#' @param useCoverage if TRUE then plot the coverage strand information, 
#' otherwise plot the number of reads strand information. FALSE by default
#' @param ... used to pass parameters to facet_wrap during plotting
#' 
#' @return The plot will be returned as a standard ggplot2 object
#' 
#' @seealso \code{\link{getWinFromBamFile}},  \code{\link{plotHist}}
#'
#' @importFrom dplyr select mutate distinct one_of starts_with
#' @importFrom stats pnorm
#' @importFrom stringr str_extract
#' @importFrom ggplot2 ggplot geom_point aes_string labs 
#' @importFrom ggplot2 theme_bw theme facet_wrap geom_line ggsave
#' @importFrom grid unit
#' @examples
#' \dontrun{
#' bamfilein = system.file("extdata","s1.sorted.bam",package = "strandCheckR")
#' windows <- getWinFromBamFile(file = bamfilein)
#' plotWin(windows)
#' 
#' # Change point colour using ggplot2
#' library(ggplot2)
#' plotWin(bamWindows) + 
#'   scale_colour_manual(values = rgb(seq(0, 1, length.out = 4), 0, 0))}
#' 
#' @export
#'
plotWin <- function(windows, split=c(10,100,1000), threshold=c(0.6,0.7,0.8,0.9),
                    save=FALSE, file="win.pdf", facets = NULL, 
                    useCoverage=FALSE, ...){

    # The initial checks for appropriate input
    reqWinCols <- c("NbPositive", "NbNegative", "CovPositive", "CovNegative", 
                    "MaxCoverage")
    stopifnot(all(reqWinCols %in% colnames(windows)))
    stopifnot(is.numeric(split))
    stopifnot(is.logical(c(save, useCoverage)))
    allows_facet_wrap <- setdiff(colnames(windows),c(reqWinCols,"Start"))
    facets <- intersect(facets,allows_facet_wrap) 

    # Check to see if we have SE or PE
    readType <- ifelse("Type" %in% colnames(windows), "PE", "SE")
    if (readType == "SE") windows$Type <- "R1"

    # Calculate the proportion of reads for the + strand, based on either 
    # coverage or the number of reads
    keepCols <- c("Nb", "Cov")[useCoverage + 1]
    windows <- as.data.frame(windows)
    windows <- dplyr::select(windows, one_of(c("MaxCoverage", facets)), 
                            starts_with(keepCols))
    names(windows) <- str_extract(names(windows),
                                paste0(c("MaxCoverage",facets,"Pos","Neg"),
                                        collapse = "|"))
    windows$NbReads <- windows$Pos + windows$Neg
    windows$PositiveProportion <- windows$Pos / windows$NbReads

    # Annotate windows by the level of coverage in the each window
    if (length(split)==0){
        windows$group <- as.factor("all")
    } 
    else{
        covBreaks <- unique(c(0, split, max(windows$MaxCoverage)))
        nBreaks <- length(covBreaks)
        breakMat <- cbind(covBreaks[-nBreaks], covBreaks[-1])
        covLabels <- apply(breakMat, MARGIN = 1, FUN = function(x){
            paste0("(", x[1], ",", x[2], "]")})
        covLabels[nBreaks - 1] <- gsub("\\(([0-9]+),.+", "> \\1", 
                                    covLabels[nBreaks-1])
        windows$group <- cut(windows$MaxCoverage, breaks = covBreaks, 
                            include.lowest = TRUE, labels = covLabels)
    }

    # Remove duplicated points for faster plotting
    windows <- mutate(windows, NbReads = round(NbReads, -1), 
                    PositiveProportion = round(PositiveProportion, 2))
    windows <- distinct(dplyr::select(windows, 
                c("NbReads", "PositiveProportion", "group", facets)))


    # Generating the threshold lines
    maxReads <- max(windows$NbReads)
    ThresholdP <- data.frame("NbReads" = c(), "PositiveProportion" = c(), 
                            "Threshold"= c())
    ThresholdN <- data.frame("NbReads" = c(), "PositiveProportion" = c(), 
                            "Threshold"= c())
    nbSampling <- 1000
    x <- floor(seq(1, sqrt(maxReads), length.out = nbSampling)^2)
    for (t in threshold){
        tP = log(t/(1-t))
        positiveReadsT <- sapply(x,function(N){
            p = seq(round(N*t),N,1) # Number of positive reads
            pP = p/N # Proportion of positive reads
            mP = log(pP/(1-pP)) # Mean prop-pos-reads (logit scale)
            sdP = sqrt(1/(N*pP*(1-pP))) # SD prop-pos-reads
            #pBinom <- pbinom(t*N,size = N, prob = pP)
            pNorm <- pnorm(tP,mean = mP, sd = sdP)
            aNorm <- which(pNorm <= 0.05)[1]
            return(p[aNorm])
        })
        tP <- data.frame("NbReads" = x, 
                        "PositiveProportion" = positiveReadsT/(x), 
                        "Threshold"= paste0(t)) 
        ThresholdP <- rbind(ThresholdP,tP)
        tN <- data.frame("NbReads" = x, 
                        "PositiveProportion" = 1-positiveReadsT/(x), 
                        "Threshold"= paste0(t)) 
        ThresholdN <- rbind(ThresholdN,tN)
    }

    # Make the basic plot
    xlab <- c("Number of reads", "Number of aligned bases")[useCoverage + 1]
    g <- ggplot() +
    geom_point(data = windows, aes_string(x = "NbReads", 
                                        y = "PositiveProportion", 
                                        colour = "group")) +
    geom_line(data = ThresholdP, aes_string(x = "NbReads", 
                                            y = "PositiveProportion", 
                                            linetype = "Threshold")) +
    geom_line(data = ThresholdN, aes_string(x = "NbReads", 
                                            y = "PositiveProportion", 
                                            linetype = "Threshold")) +
    labs(x = xlab,
        y = "Proportion of Reads on '+' Strand", 
        colour = "Max Coverage") + 
    theme_bw() +
    theme(plot.margin = unit(c(0.02,0.04,0.03,0.02), "npc"))

    if (length(facets) > 0) {
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
        dotArgs <- list(...)
        allowed <- names(formals(ggsave))
        keepArgs <- names(dotArgs) %in% setdiff(allowed, 
                                                c("filename","plot","..."))
        argList <- c(list(filename = file, plot = g), dotArgs[keepArgs])
        do.call(ggsave,argList)
    }
    return(g)
}
