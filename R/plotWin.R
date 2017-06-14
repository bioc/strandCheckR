#' @title Plot the number of reads vs positive proportion of the input windows.
#'
#' @description Plot the number of reads vs positive proportion of the input windows.

#' @param windows data frame containing the positive proportion of each window, the number of reads of the window and the maximum coverage of that window.
#' Windows can be get by calling the function getWin.
#'
#' @param group an integer vector that specifies how you want to partition the windows based on the maximum coverage. By default group = c(10,100,1000), which means that your windows will be parition into 4 groups, those have maximum coverage < 10, from 10 to 100, from 100 to 1000, and > 1000
#'
#' @param threshold a real vector that specifies which threshold lines you want to draw on the plot. The positive window above the threshold line (or negative window below the threshold line) will be kept if we use that threshold to filter the data. By default, the thresholds 0.6, 0.7, 0.8, 0.9 will be plotted.
#'
#' @param pvalue 0.05 by default, any window has pvalue below it in the test comparing with a threshold will be kept
#' @param save if TRUE, then the plot will be save into the file given by file parameter
#' @param file the file name to save to plot
#' @seealso getWin, getWinPairs, plotHist
#'
#' @examples
#'
#' #for single end bam file
#' windows <- getWin(bamfilein = "data/s1.chr1.bam")
#' plotWin(windows)
#' #for paired end bamfile
#' windowsP <- getWinPairs(bamfilein = "data/120.10.bam")
#' plotWin(windowsP)
#'
#' @export
#'

plotWin <- function(windows,threshold=c(0.6,0.7,0.8,0.9),pvalue=0.05,save=FALSE,file="win.pdf",facet_wrap_chromosomes=FALSE){

  if ("Type" %in% colnames(windows)){
    if (facet_wrap_chromosomes==TRUE){
      windowsReduced <- windows %>% dplyr::mutate(NbPositiveReads = round(NbPositiveReads, -1), NbNegativeReads = round(NbNegativeReads, -1)) %>%
        dplyr::distinct(NbPositiveReads, NbNegativeReads, Type, Chr)
      windowsReduced <- dplyr::mutate(windowsReduced,"NbReads" = NbPositiveReads + NbNegativeReads,"PositiveProportion" = NbPositiveReads/(NbPositiveReads+NbNegativeReads)) %>% dplyr::select(c(Type,Chr,NbReads,PositiveProportion))
    }
    else{
      windowsReduced <- windows %>% dplyr::mutate(NbPositiveReads = round(NbPositiveReads, -1), NbNegativeReads = round(NbNegativeReads, -1)) %>%
        dplyr::distinct(NbPositiveReads, NbNegativeReads, Type)
      windowsReduced <- dplyr::mutate(windowsReduced,"NbReads" = NbPositiveReads + NbNegativeReads,"PositiveProportion" = NbPositiveReads/(NbPositiveReads+NbNegativeReads)) %>% dplyr::select(c(Type,NbReads,PositiveProportion))
    }
  }
  else{
    if (facet_wrap_chromosomes==TRUE){
      windowsReduced <- windows %>% dplyr::mutate(NbPositiveReads = round(NbPositiveReads, -1), NbNegativeReads = round(NbNegativeReads, -1)) %>%
        dplyr::distinct(NbPositiveReads, NbNegativeReads,Chr)
      windowsReduced <- dplyr::mutate(windowsReduced,"NbReads" = NbPositiveReads + NbNegativeReads,"PositiveProportion" = NbPositiveReads/(NbPositiveReads+NbNegativeReads)) %>% dplyr::select(c(Chr,NbReads,PositiveProportion))
      }
    else{
      windowsReduced <- windows %>% dplyr::mutate(NbPositiveReads = round(NbPositiveReads, -1), NbNegativeReads = round(NbNegativeReads, -1)) %>%
        dplyr::distinct(NbPositiveReads, NbNegativeReads)
      windowsReduced <- dplyr::mutate(windowsReduced,"NbReads" = NbPositiveReads + NbNegativeReads,"PositiveProportion" = NbPositiveReads/(NbPositiveReads+NbNegativeReads)) %>% dplyr::select(c(NbReads,PositiveProportion))
    }
  }
  rm(windows)

  ThresholdP <- data.frame("NbReads" = c(), "PositiveProportion" = c(), "Threshold"= c())
  ThresholdN <- data.frame("NbReads" = c(), "PositiveProportion" = c(), "Threshold"= c())
  for (t in threshold){
    positiveReadsT <- sapply(1:10000,function(N){
      tP = log(t/(1-t))
      p = seq(round(N*t),N,1)
      pP = p/N
      mP = log(pP/(1-pP))
      sdP = sqrt(1/(N*pP*(1-pP)))

      pNorm <- pnorm(tP,mean = mP, sd = sdP)
      #pBinom <- pbinom(0.9*N,size = N, prob = pP)
      aNorm <- which(pNorm <= 0.05)[1]
      return(p[aNorm])
    })
    tP <- data.frame("NbReads" = 1:10000, "PositiveProportion" = positiveReadsT/(1:10000), "Threshold"= paste0(t)) #%>% rbind(data.frame("NbPositiveReads"=max(windowsReduced$NbPositiveReads),"NbNegativeReads"=round(max(windowsReduced $NbPositiveReads)/t)-max(windowsReduced $NbPositiveReads),"Threshold"=paste0(t)))
    ThresholdP <- rbind(ThresholdP,tP)
    tN <- data.frame("NbReads" = 1:10000, "PositiveProportion" = 1-positiveReadsT/(1:10000), "Threshold"= paste0(t)) #%>% rbind(data.frame("NbPositiveReads"=round(max(windowsReduced $NbNegativeReads)/t)-max(windowsReduced $NbNegativeReads),"NbNegativeReads"=max(windowsReduced $NbNegativeReads),"Threshold"=paste0(t)))
    ThresholdN <- rbind(ThresholdN,tN)
  }
  ThresholdP <- rbind(ThresholdP,data.frame("NbReads" = max(windowsReduced$NbReads), "PositiveProportion" = threshold, "Threshold"=paste0(threshold)))
  ThresholdN <- rbind(ThresholdN,data.frame("NbReads" = max(windowsReduced$NbReads), "PositiveProportion" = (1-threshold), "Threshold"=paste0(threshold)))
  if ("Type" %in% colnames(windows)){
    gg <- ggplot2::ggplot()
    if (facet_wrap_chromosomes==TRUE){
      gg <- gg + ggplot2::geom_point(data = windowsReduced, ggplot2::aes(x = NbReads, y = PositiveProportion)) +
        ggplot2::facet_wrap(~Type+Chr)
    }
    else{
      gg <- gg + ggplot2::geom_point(data = windowsReduced, ggplot2::aes(x = NbReads, y = PositiveProportion)) +
        ggplot2::facet_wrap(~Type)
    }
    gg <- gg + ggplot2::geom_line(data = ThresholdP, ggplot2::aes(x = NbReads, y = PositiveProportion, colour = Threshold))
    gg + ggplot2::geom_line(data = ThresholdN, ggplot2::aes(x = NbReads, y = PositiveProportion, colour = Threshold)) + ggplot2::theme_bw()
    if (save==TRUE) {
      message("The plot will be saved to the file ",file)
      ggplot2::ggsave(filename = file)
    }
    else{
      gg
    }
  }
  else{
    gg <- ggplot2::ggplot()
    if (facet_wrap_chromosomes==TRUE){
      gg <- gg + ggplot2::geom_point(data = windowsReduced, ggplot2::aes(x = NbReads, y = PositiveProportion))+ ggplot2::facet_wrap(~Chr)
    }
    else{
      gg <- gg + ggplot2::geom_point(data = windowsReduced, ggplot2::aes(x = NbReads, y = PositiveProportion))
    }
    gg <- gg + ggplot2::geom_line(data = ThresholdP, ggplot2::aes(x = NbReads, y = PositiveProportion, colour = Threshold))
    gg <- gg + ggplot2::geom_line(data = ThresholdN, ggplot2::aes(x = NbReads, y = PositiveProportion, colour = Threshold)) + ggplot2::theme_bw()
    if (save==TRUE){
      message("The plot will be saved to the file ",file)
      ggplot2::ggsave(filename = file, plot = gg)
    }
    else{
     gg
    }
  }
}
