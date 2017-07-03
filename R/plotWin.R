#' @title Plot the number of reads vs positive proportion of the input windows data frame.
#'
#' @description Plot the number of reads vs positive proportion of the input windows data frame.

#' @param windows data frame containing the number of positive/negative reads for each window.
#' Windows can be get by the function \code{getWinFromBamFile} (for single end bam file) or \code{getWinFromPairedBamFile} (for paired end bam file).
#'
#' @param group an integer vector that specifies how you want to partition the windows based on the number of reads. By default \code{group} = c(10,100,1000), which means that your windows will be parition into 4 groups, those have number of reads < 10, from 10 to 100, from 100 to 1000, and > 1000
#'
#' @param threshold a real vector that specifies which threshold lines you want to draw on the plot. The positive window above the threshold line (or negative window below the threshold line) will be kept if we use that threshold to filter the data. By default, the thresholds 0.6, 0.7, 0.8, 0.9 will be plotted.
#'
#' @param pvalue 0.05 by default, any window has pvalue below it in the test comparing with a threshold will be kept.
#' @param save if TRUE, then the plot will be save into the file given by \code{file} parameter
#' @param file the file name to save to plot
#' @param facet_wrap_chromosomes if TRUE, then the plots will be splitted by chromosomes. FALSE by default
#' @seealso getWinFromBamFile, getWinFromPairedBamFile, plotHist
#'
#' @examples
#'
#' #for single end bam file
#' bamfilein = system.file("data/s1.chr1.bam",package = "rnaCleanR")
#' windows <- getWinFromBamFile(bamfilein)
#' plotWin(windows)
#' #for paired end bamfile
#' bamfilepair = system.file("data/120.10.bam",package = "rnaCleanR")
#' windowsP <- getWinFromPairedBamFile(bamfilein = "data/120.10.bam")
#' plotWin(windowsP)
#'
#' @export
#'

plotWin <- function(windows,group=c(10,100,1000),threshold=c(0.6,0.7,0.8,0.9),pvalue=0.05,save=FALSE,file="win.pdf",facet_wrap_chromosomes=FALSE){
  coverage <- "MaxCoverage" %in% colnames(windows)
  if (!facet_wrap_chromosomes) windows <- dplyr::select(windows,-Chr)
  windows <- dplyr::select(windows,-Start) 
  windows <- dplyr::mutate(windows,"NbReads"=NbPositiveReads+NbNegativeReads,"PositiveProportion" = NbPositiveReads/(NbPositiveReads+NbNegativeReads)) %>%
      dplyr::select(-c(NbPositiveReads,NbNegativeReads))
  
  if (coverage){
    if (length(group)==0){
      leg <- "all"
    } else{
      leg <- paste0("<",group[1])
      for (i in seq_along(group[-1])){
        leg <- c(leg,paste0(group[i],"-",group[i+1]))
      }
      leg <- c(leg,paste0(">",group[length(group)]))
    }
    x <- lapply(seq_along(group),function(i){which(group[i]<windows$MaxCoverage)})
    G <- rep(leg[1],nrow(windows))
    for (i in seq_along(group)){
      G[x[[i]]] <- leg[i+1]
    }
    windows$MaxCoverage = G  
  }
  
  if ("Type" %in% colnames(windows)){
    if (facet_wrap_chromosomes==TRUE){
      if (coverage){
        windowsReduced <- windows %>% dplyr::mutate(NbReads = round(NbReads, -1), PositiveProportion = round(PositiveProportion, 3)) %>%
          dplyr::distinct(Chr,Type,NbReads,MaxCoverage,PositiveProportion) 
      } else{
        windowsReduced <- windows %>% dplyr::mutate(NbReads = round(NbReads, -1), PositiveProportion = round(PositiveProportion, 3)) %>%
          dplyr::distinct(Chr,Type,NbReads, PositiveProportion)  
      }
    } else{
      if (coverage){
        windowsReduced <- windows %>% dplyr::mutate(NbReads = round(NbReads, -1), PositiveProportion = round(PositiveProportion, 3)) %>%
          dplyr::distinct(Type,NbReads,MaxCoverage,PositiveProportion)
      } else{
        windowsReduced <- windows %>% dplyr::mutate(NbReads = round(NbReads, -1), PositiveProportion = round(PositiveProportion, 3)) %>%
          dplyr::distinct(Type,NbReads, PositiveProportion)  
      }
    }
  } else{
    if (facet_wrap_chromosomes==TRUE){
      if (coverage){
        windowsReduced <- windows %>% dplyr::mutate(NbReads = round(NbReads, -1), PositiveProportion = round(PositiveProportion, 3)) %>%
          dplyr::distinct(Chr,NbReads,MaxCoverage,PositiveProportion)  
      } else{
        windowsReduced <- windows %>% dplyr::mutate(NbReads = round(NbReads, -1), PositiveProportion = round(PositiveProportion, 3)) %>%
          dplyr::distinct(Chr,NbReads, PositiveProportion)
      }
    } else{
      if (coverage){
        windowsReduced <- windows %>% dplyr::mutate(NbReads = round(NbReads, -1), PositiveProportion = round(PositiveProportion, 3)) %>%
          dplyr::distinct(NbReads,MaxCoverage,PositiveProportion)  
      } else{
        windowsReduced <- windows %>% dplyr::mutate(NbReads = round(NbReads, -1), PositiveProportion = round(PositiveProportion, 3)) %>%
          dplyr::distinct(NbReads, PositiveProportion)
      }
      
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
  if ("Type" %in% colnames(windowsReduced)){
    gg <- ggplot2::ggplot()
    if (coverage){
      gg <- gg + ggplot2::geom_point(data = windowsReduced, ggplot2::aes(x = NbReads, y = PositiveProportion, colour = MaxCoverage))  
      gg <- gg + ggplot2::geom_line(data = ThresholdP, ggplot2::aes(x = NbReads, y = PositiveProportion, linetype = Threshold))
      gg <- gg + ggplot2::geom_line(data = ThresholdN, ggplot2::aes(x = NbReads, y = PositiveProportion, linetype = Threshold))
    } else{
      gg <- gg + ggplot2::geom_point(data = windowsReduced, ggplot2::aes(x = NbReads, y = PositiveProportion))
      gg <- gg + ggplot2::geom_line(data = ThresholdP, ggplot2::aes(x = NbReads, y = PositiveProportion, colour = Threshold))
      gg <- gg + ggplot2::geom_line(data = ThresholdN, ggplot2::aes(x = NbReads, y = PositiveProportion, colour = Threshold))
    }
    if (facet_wrap_chromosomes==TRUE){
        gg <- gg + ggplot2::facet_wrap(~Type+Chr)  
    } else{
        gg <- gg +ggplot2::facet_wrap(~Type)  
    }
    gg <- gg + ggplot2::theme_bw()
    if (coverage){
      gg <- gg + ggplot2::labs(x = "Number of Bases", y = "Positive Proportion")
    } else{
      gg <- gg + ggplot2::labs(x = "Number of Reads", y = "Positive Proportion")
    }
    if (save==TRUE) {
      message("The plot will be saved to the file ",file)
      ggplot2::ggsave(filename = file)
    } else{
      gg
    }
  } else{
    gg <- ggplot2::ggplot()
    if (coverage){
      gg <- gg + ggplot2::geom_point(data = windowsReduced, ggplot2::aes(x = NbReads, y = PositiveProportion, colour = MaxCoverage))  
      gg <- gg + ggplot2::geom_line(data = ThresholdP, ggplot2::aes(x = NbReads, y = PositiveProportion, linetype = Threshold))
      gg <- gg + ggplot2::geom_line(data = ThresholdN, ggplot2::aes(x = NbReads, y = PositiveProportion, linetype = Threshold))
    } else{
      gg <- gg + ggplot2::geom_point(data = windowsReduced, ggplot2::aes(x = NbReads, y = PositiveProportion))
      gg <- gg + ggplot2::geom_line(data = ThresholdP, ggplot2::aes(x = NbReads, y = PositiveProportion, colour = Threshold))
      gg <- gg + ggplot2::geom_line(data = ThresholdN, ggplot2::aes(x = NbReads, y = PositiveProportion, colour = Threshold))
    }
    if (facet_wrap_chromosomes==TRUE){
      gg <- gg +ggplot2::facet_wrap(~Chr)  
    }
    gg <- gg + ggplot2::theme_bw()
    if (coverage){
      gg <- gg + ggplot2::labs(x = "Number of Bases", y = "Positive Proportion")
    } else{
      gg <- gg + ggplot2::labs(x = "Number of Reads", y = "Positive Proportion")
    }
    if (save==TRUE) {
      message("The plot will be saved to the file ",file)
      ggplot2::ggsave(filename = file)
    } else{
      gg
    }
  }
}
