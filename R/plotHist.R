#' @title Plot the histogram of positive proportions 
#'
#' @description Plot the histogram of positive proportions of the input 
#' histogram data frame
#'
#' @param windows data frame containing the strand information of the sliding 
#' windows. Windows can be obtained using the function \code{getWinFromBamFile}.
#' @param save if TRUE, then the plot will be save into the file given by 
#' \code{file} parameter
#' @param file the file name to save to plot
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
#' @param heatmap if TRUE, then use heat map to plot the histogram, otherwise 
#' use barplot. FALSE by default.
#' @param ... used to pass parameters to facet_wrap
#' 
#' @return If \code{heatmap=FALSE}: a ggplot object
#' 
#' @seealso \code{\link{getWinFromBamFile}}, \code{\link{plotWin}}
#'
#' @examples
#' bamfilein = system.file("extdata","s1.sorted.bam",package = "strandCheckR")
#' win  <- getWinFromBamFile(file = bamfilein)
#' plotHist(win)
#' 
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_bar geom_tile
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw theme element_blank
#' @importFrom grid unit
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 ggsave ggtitle scale_fill_gradient
#' @importFrom magrittr set_colnames
#' @importFrom dplyr mutate select one_of starts_with bind_rows filter
#' @importFrom stringr str_extract
#' @importFrom reshape2 dcast melt
#' @export
plotHist <- function(windows, save=FALSE, file = "hist.pdf", group_by = NULL, 
                    normalize_by = NULL, split=c(10,100,1000), breaks = 100,
                    useCoverage=FALSE, heatmap = FALSE, ...){
    histWin <- summarizeHist(windows, split = split, breaks = breaks, 
                            useCoverage = useCoverage, group_by = group_by,
                            normalize_by = normalize_by)
    # The initial checks for appropriate input
    reqWinCols <- c("PosStrandProp","ReadCountProp")
    stopifnot(all(reqWinCols %in% colnames(histWin)))
    stopifnot(is.logical(save))
    allows_facet_wrap <- setdiff(colnames(histWin),c(reqWinCols,"Coverage"))
    group_by <- intersect(group_by,allows_facet_wrap) 

    # Make the plot
    if (heatmap==FALSE){
        g <- ggplot(histWin, aes_string(x = "PosStrandProp", 
                                        y = "ReadCountProp", 
                                        fill = "Coverage")) + 
            geom_bar(stat = "identity") + 
            labs(x = "Proportion of Reads on '+' Strand",
                y = "Proportion of Windows") +
            theme_bw() +
            theme(plot.margin = unit(c(0.02,0.04,0.03,0.02), "npc"))  
        if (length(group_by)>0) {

            # Get any facet arguments from dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(facet_wrap))
            keepArgs <- names(dotArgs) %in% setdiff(allowed, "facets")
            argList <- c(list(facets = group_by), dotArgs[keepArgs])
            myFacets <- do.call(facet_wrap, argList)
            g <- g + myFacets
        }
    } else{
        cov <- unique(histWin$Coverage)
        g <- list()
        if (!("File" %in% colnames(histWin))){
            histWin$File <- "File"
        }
        group_by <- group_by[group_by!="File"]
        if (length(group_by)>0) {
            # Get any facet arguments from dotArgs that have been set manually
            dotArgs <- list(...)
            allowed <- names(formals(facet_wrap))
            keepArgs <- names(dotArgs) %in% setdiff(allowed, "facets")
            argList <- c(list(facets = group_by), dotArgs[keepArgs])
            myFacets <- do.call(facet_wrap, argList)
        }
        for (i in seq_along(cov)){
            l <- filter(histWin,Coverage==cov[i]) 
            g[[i]] <- ggplot(l, aes_string(x="PosStrandProp", y="File",
                                        fill="ReadCountProp")) + 
                geom_tile() +
                theme_bw() +
                theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) +
                scale_fill_gradient(low="white",high="red",na.value = "white")
            if (length(group_by)>0){
                g[[i]] <- g[[i]] + myFacets
            }
            if (length(cov)>1) {
                g[[i]] <- g[[i]] + ggtitle(paste0("Coverage ",cov[i]))
            }
        }
        g <- grid.arrange(grobs=g, nrow = length(cov))
    }

    if (save==TRUE){
        message("The plot will be saved to the file ",file)
        dotArgs <- list(...)
        allowed <- names(formals(ggsave))
        keepArgs <- names(dotArgs) %in% setdiff(allowed, c("filename","plot"))
        argList <- c(list(filename = file, plot = g), dotArgs[keepArgs])
        do.call(ggsave,argList)
    }
    if (heatmap==FALSE) return(g)
}

