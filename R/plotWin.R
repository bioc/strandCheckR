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
#' @seealso getWin, plotHist
#' 
#' @examples 
#' 
#' #for single end bam file
#' windows <- getWin(bamfilein = "data/s1.chr1.bam",readLength=50)
#' plotWin(windows)
#' #for paired end bamfile
#' windowsP <- getWinPairs(bamfilein = "data/120.10.bam",readLength=100)
#' plotWin(windowsP)
#' 
#' @export
#'

plotWin <- function(windows,group=c(10,100,1000),threshold=c(0.6,0.7,0.8,0.9),pvalue=0.05){
  if (length(group)==0){
    leg <- "all"
  }
  else{
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
  
  if ("Type" %in% colnames(windows)){
    windowsReduced <- windows %>% dplyr::mutate(NbPositiveReads = round(NbPositiveReads, -1), NbNegativeReads = round(NbNegativeReads, -1)) %>% 
      dplyr::distinct(NbPositiveReads, NbNegativeReads, MaxCoverage,Type)  
  }
  else{
    windowsReduced <- windows %>% dplyr::mutate(NbPositiveReads = round(NbPositiveReads, -1), NbNegativeReads = round(NbNegativeReads, -1)) %>% 
      dplyr::distinct(NbPositiveReads, NbNegativeReads, MaxCoverage)  
  }
  windowsReduced <- dplyr::mutate(windowsReduced,"NbReads" = NbPositiveReads + NbNegativeReads,"PositiveProportion" = NbPositiveReads/(NbPositiveReads+NbNegativeReads))
  # for (t in threshold){
  #   pP = seq(0.5,1,0.005)
  #   mP = log(pP/(1-pP))
  #   tP = log(t/(1-t))
  #   pN = seq(0,0.5,0.005)
  #   mN = log(pN/(1-pN))
  #   tN = log((1-t)/t)
  #   NbReads = 50000
  #   propor <- sapply(1:NbReads,function(N){
  #     sdP = sqrt(1/(N*pP*(1-pP)))
  #     xP <- pnorm(tP,mean=mP,sd=sdP)
  #     aP <- which(xP<=0.05)[1]
  #     return(pP[aP])
  #   }) 
  #   #p <- seq(t,1,0.0005)
  #   #n = ((qnorm(1-pvalue)/sqrt(p*(1-p))/(log(p/(1-p))-log(t/(1-t))))^2)
  #   thresholdP <- rbind(thresholdP,data.frame("Proportion" = c(propor,t), "NbReads" = c(1:NbReads,Inf), "Threshold" = paste0(t)))
  # }
  # thresholdN <- thresholdP
  # thresholdN$Proportion <- 1-thresholdN$Proportion
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
    gg <- gg + ggplot2::geom_point(data = windowsReduced, ggplot2::aes(x = NbReads, y = PositiveProportion, colour = MaxCoverage)) 
    gg <- gg + ggplot2::geom_line(data = ThresholdP, ggplot2::aes(x = NbReads, y = PositiveProportion, linetype = Threshold))
    gg + ggplot2::geom_line(data = ThresholdN, ggplot2::aes(x = NbReads, y = PositiveProportion, linetype = Threshold)) + ggplot2::theme_bw() + ggplot2::facet_wrap(~Type)
  }
  else{
    gg <- ggplot2::ggplot()
    gg <- gg + ggplot2::geom_point(data = windowsReduced, ggplot2::aes(x = NbReads, y = PositiveProportion, colour = MaxCoverage)) 
    gg <- gg + ggplot2::geom_line(data = ThresholdP, ggplot2::aes(x = NbReads, y = PositiveProportion, linetype = Threshold))
    gg + ggplot2::geom_line(data = ThresholdN, ggplot2::aes(x = NbReads, y = PositiveProportion, linetype = Threshold)) + ggplot2::theme_bw()
  }
}