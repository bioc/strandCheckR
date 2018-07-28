#' @title Summarize the histogram of strand proportions from the input windows 
#' data frame
#'
#' @description Summarize the histogram of positive proportions from the input 
#' windows obtained from the function \code{getWinFromBamFile}
#'
#' @param windows data frame containing the strand information of the sliding 
#' windows. Windows can be obtained using the function \code{getWinFromBamFile}.
#' @param split an integer vector that specifies how you want to partition the 
#' windows based on the coverage. By default \code{split} = c(10,100,1000), 
#' which means that your windows will be partitionned into 4 groups, those have 
#' coverage < 10, from 10 to 100, from 100 to 1000, and > 1000
#' @param breaks an integer giving the number of bins for the histogram
#' @param useCoverage if TRUE then plot the coverage strand information, 
#' otherwise plot the number of reads strand information. FALSE by default
#' @param group_by the column names of windows that will be used to group 
#' the data
#' @param normalize_by the column names of windows that will be used to 
#' normalize the read count or read coverage into proportion
#' 
#' @return a dataframe object
#' 
#' @seealso getWinFromBamFile, plotHist, plotWin
#' 
#' @importFrom magrittr set_colnames
#' @importFrom dplyr mutate select one_of starts_with bind_rows
#' @importFrom stringr str_extract
#' @importFrom reshape2 dcast melt

summarizeHist <- function(
    windows, split = c(10, 100, 1000), breaks = 100, useCoverage = FALSE, 
    group_by = NULL, normalize_by = NULL
    ) 
{   
    # The initial checks for appropriate input
    reqWinCols <- c(
        "NbPos", "NbNeg", "CovPos", "CovNeg", "MaxCoverage"
        )
    stopifnot(all(reqWinCols %in% colnames(windows)))
    stopifnot(length(split) == 0 || is.numeric(split))
    stopifnot(is.logical(useCoverage))
    
    # Check to see if we have SE or PE
    readType <- ifelse("Type" %in% colnames(windows), "PE", "SE")
    if (readType == "SE") windows$Type <- "R1"
    
    allows_group_by <- setdiff(colnames(windows), c(reqWinCols, "Start", "End"))
    group_by <- intersect(group_by, allows_group_by)
    if (length(group_by) > 0){
        message("Windows are grouped by",group_by,"\n")
    }
    # Calculate the proportion of reads for the + strand, based on either 
    # coverage or the number of reads
    keepCols <- c("Nb", "Cov")[useCoverage + 1]
    windows <- as.data.frame(windows)
    windows <- dplyr::select(
        windows, one_of(c("MaxCoverage", group_by)), starts_with(keepCols)
        )
    names(windows) <- str_extract(
        names(windows), 
        paste0(c("MaxCoverage", "Pos", "Neg", group_by), collapse = "|")
        )
    if (useCoverage == FALSE) windows <- filter(windows, Pos + Neg > 0)
    windows$PosProp <- windows$Pos/(windows$Pos + windows$Neg)
    
    # Annotate windows by the level of coverage in the each window
    if (length(split) == 0) {
        windows$group <- as.factor("all")
    } else {
        covBreaks <- unique(c(0, split, max(windows$MaxCoverage)))
        nBreaks <- length(covBreaks)
        breakMat <- cbind(covBreaks[-nBreaks], covBreaks[-1])
        covLabels <- apply(
            breakMat, MARGIN = 1, 
            FUN = function(x) {paste0("(", x[1], ",", x[2], "]")}
            )
        covLabels[nBreaks - 1] <- 
            gsub("\\(([0-9]+),.+", "> \\1", covLabels[nBreaks - 1])
        windows$group <- cut(
            windows$MaxCoverage, breaks = covBreaks, include.lowest = TRUE, 
            labels = covLabels
            )
    }
    
    # Group the stranded proportions into bins
    windows$PosProp <- cut(
        windows$PosProp, breaks = seq(0, 1, length.out = breaks + 2), 
        include.lowest = TRUE
        )
    # Tally the bins by coverage
    formula <- paste0(
        paste0(c(group_by, "PosProp"), collapse = "+"), " ~ group"
        )
    windows <- dcast(
        windows, formula, fun.aggregate = length, value.var = "PosProp"
        )
    # Melt for easy plotting with ggplot2
    windows <- melt(
        windows, id.vars = c(group_by, "PosProp"), variable.name = "Coverage", 
        value.name = "ReadCountProp"
        )
    windows$PosProp <- (as.integer(windows$PosProp) - 1)/breaks
    # Convert the numbers of windows into proportions keeping R1 & R2 separate
    normalize_by <- intersect(normalize_by, group_by)
    if (length(normalize_by) > 0) {
        message("Windows are normalized by", normalize_by,"\n")
        windows <- lapply(
            split(windows, f = windows[, normalize_by]), 
            function(x) {
                x$ReadCountProp <- x$ReadCountProp/sum(x$ReadCountProp)
                x
                }
            )
        windows <- bind_rows(windows)
    } else {
        windows$ReadCountProp <- 
            windows$ReadCountProp/sum(windows$ReadCountProp)
    }
    return(windows)
}

