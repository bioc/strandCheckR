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
#' windows <- data.frame("Proportion"=runif(1000, min=0, max=1),"NbReads" = sample(1:5000,1000,replace=TRUE)) %>% mutate("MaxCoverage"=sample(NbReads, 1000, replace=TRUE))
#' windows <- windows[order(windows$NbReads),]
#' plotWin(windows,group=c(10,100,1000))
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
  
  thresholdP <- data.frame("Proportion" = c(), "NbReads" = c(), "Threshold"=c())
  
  for (t in threshold){
    pP = seq(0.5,1,0.005)
    mP = log(pP/(1-pP))
    tP = log(t/(1-t))
    pN = seq(0,0.5,0.005)
    mN = log(pN/(1-pN))
    tN = log((1-t)/t)
    NbReads = 50000
    propor <- sapply(1:NbReads,function(N){
      sdP = sqrt(1/(N*pP*(1-pP)))
      xP <- pnorm(tP,mean=mP,sd=sdP)
      aP <- which(xP<=0.05)[1]
      return(pP[aP])
    }) 
    #p <- seq(t,1,0.0005)
    #n = ((qnorm(1-pvalue)/sqrt(p*(1-p))/(log(p/(1-p))-log(t/(1-t))))^2)
    thresholdP <- rbind(thresholdP,data.frame("Proportion" = c(propor,t), "NbReads" = c(1:NbReads,Inf), "Threshold" = paste0(t)))
  }
  thresholdN <- thresholdP
  thresholdN$Proportion <- 1-thresholdN$Proportion
  
  if ("Type" %in% colnames(windows)){
    windowsReduced <- windows %>% mutate(NbReads = round(NbReads, -1), Proportion = round(Proportion, 2)) %>% distinct(Proportion, NbReads, MaxCoverage,Type)  
    gg <- ggplot()
    gg <- gg + geom_point(data = windowsReduced, aes(x = NbReads, y = Proportion, colour = MaxCoverage)) + facet_wrap(~Type)
    gg <- gg + geom_line(data = thresholdP, aes(x = NbReads, y = Proportion, linetype = Threshold)) 
    gg + geom_line(data = thresholdN, aes(x = NbReads, y = Proportion, linetype = Threshold)) + xlim(min(windowsReduced$NbReads),max(windowsReduced$NbReads)) + ylim(0,1) + theme_bw()
  }
  else{
    windowsReduced <- windows %>% mutate(NbReads = round(NbReads, -1), Proportion = round(Proportion, 2)) %>% distinct(Proportion, NbReads, MaxCoverage)
    gg <- ggplot()
    gg <- gg + geom_point(data = windowsReduced, aes(x = NbReads, y = Proportion, colour = MaxCoverage))
    gg <- gg + geom_path(data = thresholdP, aes(x = NbReads, y = Proportion, linetype = Threshold)) 
    gg + geom_path(data = thresholdN, aes(x = NbReads, y = Proportion, linetype = Threshold)) + xlim(min(windowsReduced$NbReads),max(windowsReduced$NbReads)) + ylim(0,1) + theme_bw()
  }
  
}