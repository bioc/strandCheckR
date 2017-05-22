#' @title Plot the histogram of positive proportions of the input windows
#' 
#' @description Plot the histogram of positive proportions of the input windows
#' 
#' @param windows data frame containing the positive proportion of each window, the number of reads of the window and the maximum coverage of that window.
#' Windows can be get by calling the function getWin. 
#' 
#' @param group an integer vector that specifies how you want to partition the windows based on the maximum coverage. By default group = c(10,100,1000), which means that your windows will be parition into 4 groups, those have maximum coverage < 10, from 10 to 100, from 100 to 1000, and > 1000
#' 
#' @seealso getWin, plotWin
#' 
#' @examples
#' #for single end bam file
#' windows <- getWin(bamfilein = "data/s1.chr1.bam",readLength=50)
#' plotHist(windows)
#' #for paired end bamfile
#' windowsP <- getWinPairs(bamfilein = "data/120.10.bam",readLength=100)
#' plotHist(windowsP)
#' @export
#'
plotHist <- function(windows,group=c(10,100,1000)){
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
  breaks <- 100
  if ("Type" %in% colnames(windows)){
    histoFirst <- lapply(leg,function(l){
      a <- dplyr::filter(windows,MaxCoverage==l,Type=="First")$Proportion
      if (length(a)>0){
        hist(a,breaks=0:(breaks+1)/breaks - 1/(2*breaks),plot=FALSE)$count 
      }}) 
    null <- sapply(histoFirst,is.null)
    histoFirst <- histoFirst[!null] %>% 
      as.data.frame() %>%  
      set_colnames(leg[!null]) %>%
      dplyr::mutate("Type"="First","Proportion"=0:100) %>% 
      melt(id.vars = c("Type","Proportion")) %>% 
      set_colnames(c("Type","Proportion","MaxCoverage","Count"))
    
    histoSecond <- lapply(leg,function(l){
      a <- dplyr::filter(windows,MaxCoverage==l,Type=="Second")$Proportion
      if (length(a)>0){
        hist(a,breaks=0:(breaks+1)/breaks - 1/(2*breaks),plot=FALSE)$count  
      }}) 
    null <- sapply(histoSecond,is.null)
    histoSecond <- histoSecond[!null] %>% 
      as.data.frame() %>%  
      set_colnames(leg[!null]) %>%
      dplyr::mutate("Type"="Second","Proportion"=0:100) %>% 
      melt(id.vars = c("Type","Proportion")) %>% 
      set_colnames(c("Type","Proportion","MaxCoverage","Count"))
    
    ggplot(rbind(histoFirst,histoSecond), aes(Proportion, Count, fill=MaxCoverage, width=1))+
      facet_wrap(~Type) +
      geom_bar(stat="identity", width=.3,position="stack") + 
      ylab("Count") + 
      xlab("Positive Proportion")+ 
      theme_bw() 
  }
  else{
    histo<- lapply(leg,function(l){
      a <- dplyr::filter(windows,MaxCoverage==l)$Proportion
      if (length(a)>0){
        hist(a,breaks=0:(breaks+1)/breaks - 1/(2*breaks),plot=FALSE)$count
      }}) 
    null <- sapply(histo,is.null)
    histo <- histo[!null] %>% 
      as.data.frame() %>%  
      set_colnames(leg[!null]) %>%
      dplyr::mutate("Proportion"=0:100) %>% 
      reshape2::melt(id.vars = c("Proportion")) %>% 
      set_colnames(c("Proportion","MaxCoverage","Count"))
    ggplot(histo, aes(Proportion, Count, fill=MaxCoverage, width=1))+
      geom_bar(stat="identity", width=.3,position="stack") + 
      ylab("Count") + 
      xlab("Positive Proportion")+ 
      theme_bw() 
  }
}
