#' @title Plot the histogram of positive proportions of the input windows
#'
#' @description Plot the histogram of positive proportions of the input windows
#'
#' @param windows data frame containing the positive proportion of each window, the number of reads of the window and the maximum coverage of that window.
#' Windows can be get by calling the function getWin.
#'
#' @param group an integer vector that specifies how you want to partition the windows based on the maximum coverage. By default group = c(10,100,1000), which means that your windows will be parition into 4 groups, those have maximum coverage < 10, from 10 to 100, from 100 to 1000, and > 1000
#' @param save if TRUE, then the plot will be save into the file given by file parameter
#' @param file the file name to save to plot
#' @param facet_wrap_chromosomes if TRUE, then the plots will be splitted by chromosomes. FALSE by default
#' @seealso getWin, getWinPairs, plotWin
#'
#' @examples
#' #for single end bam file
#' windows <- getWin(bamfilein = "data/s1.chr1.bam")
#' plotHist(windows)
#' #for paired end bamfile
#' windowsP <- getWinPairs(bamfilein = "data/120.10.bam")
#' plotHist(windowsP)
#' @export
#' @importFrom magrittr set_colnames
#'
plotHist <- function(windows,group=c(10,100,1000),save=FALSE,file = "hist.pdf",facet_wrap_chromosomes=FALSE){
  windows <- dplyr::mutate(windows,"NbReads" = windows$NbPositiveReads+windows$NbNegativeReads) %>% dplyr::select(-Start)
  if (!facet_wrap_chromosomes) windows <- dplyr::select(windows,-Chr)
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
  x <- lapply(seq_along(group),function(i){which(group[i]<windows$NbReads)})
  G <- rep(leg[1],nrow(windows))
  for (i in seq_along(group)){
    G[x[[i]]] <- leg[i+1]
  }
  windows$NbReads = G
  breaks <- 100
  if ("Type" %in% colnames(windows)){
    histoFirst <- lapply(leg,function(l){
      a <- dplyr::filter(windows,NbReads==l,Type=="First") %>% dplyr::mutate("Proportion"=NbPositiveReads/(NbPositiveReads+NbNegativeReads))
      if (nrow(a)>0){
        if (facet_wrap_chromosomes){
          chromosomes <- unique(a$Chr)
          sapply(chromosomes,function(chr){
            ac <- dplyr::filter(a,Chr==chr)
            hist(ac$Proportion,breaks=0:(breaks+1)/breaks - 1/(2*breaks),plot=FALSE)$count
          }) %>% reshape2::melt()  %>% set_colnames(c("Proportion","Chr","Count")) %>% dplyr::mutate("NbReads"=l)
        }
        else{
          data.frame("Count"=hist(a$Proportion,breaks=0:(breaks+1)/breaks - 1/(2*breaks),plot=FALSE)$count,
                     "Proportion" = 0:100,
                     "NbReads" = l)
        }
      }})
    null <- sapply(histoFirst,is.null)
    histoFirst <- histoFirst[!null]
    histoFirst <- do.call(rbind,histoFirst) %>% dplyr::mutate("Type"="First")

    histoSecond <- lapply(leg,function(l){
      a <- dplyr::filter(windows,NbReads==l,Type=="Second") %>% dplyr::mutate("Proportion"=NbPositiveReads/(NbPositiveReads+NbNegativeReads))
      if (nrow(a)>0){
        if (facet_wrap_chromosomes){
          chromosomes <- unique(a$Chr)
          sapply(chromosomes,function(chr){
            ac <- dplyr::filter(a,Chr==chr)
            hist(ac$Proportion,breaks=0:(breaks+1)/breaks - 1/(2*breaks),plot=FALSE)$count
          }) %>% reshape2::melt()  %>% set_colnames(c("Proportion","Chr","Count")) %>% dplyr::mutate("NbReads"=l)
        }
        else{
          data.frame("Count"=hist(a$Proportion,breaks=0:(breaks+1)/breaks - 1/(2*breaks),plot=FALSE)$count,
                     "Proportion" = 0:100,
                     "NbReads" = l)
        }
      }})
    null <- sapply(histoSecond,is.null)
    histoSecond <- histoSecond[!null]
    histoSecond <- do.call(rbind,histoSecond) %>% dplyr::mutate("Type"="Second")
    histo <- rbind(histoFirst,histoSecond)
    if (facet_wrap_chromosomes){
      g <- ggplot2::ggplot(histo, ggplot2::aes(Proportion, Count, fill=NbReads, width=1))+
        scale_fill_discrete(breaks=leg)+
        ggplot2::geom_bar(stat="identity", width=.3,position="stack") +
        ggplot2::facet_wrap(~Type+Chr)
    }
    else{
      g <- ggplot2::ggplot(histo, ggplot2::aes(Proportion, Count, fill=NbReads, width=1))+
        scale_fill_discrete(breaks=leg)+
        ggplot2::geom_bar(stat="identity", width=.3,position="stack")+
        ggplot2::facet_wrap(~Type)
    }
     g <- g +  ggplot2::ylab("Count") +
      ggplot2::xlab("Positive Proportion")+
      ggplot2::theme_bw()
    if (save==TRUE){
      message("The plot will be saved to the file ",file)
      ggplot2::ggsave(filename = file)
    }
    else{
     g
    }
  }
  else{
    histo<- lapply(leg,function(l){
      a <- dplyr::filter(windows,NbReads==l) %>% dplyr::mutate("Proportion"=NbPositiveReads/(NbPositiveReads+NbNegativeReads))
      if (nrow(a)>0){
        if (facet_wrap_chromosomes){
          chromosomes <- unique(a$Chr)
          sapply(chromosomes,function(chr){
            ac <- dplyr::filter(a,Chr==chr)
            hist(ac$Proportion,breaks=0:(breaks+1)/breaks - 1/(2*breaks),plot=FALSE)$count
          }) %>% reshape2::melt()  %>% set_colnames(c("Proportion","Chr","Count")) %>% dplyr::mutate("NbReads"=l)
        }
        else{
          data.frame("Count"=hist(a$Proportion,breaks=0:(breaks+1)/breaks - 1/(2*breaks),plot=FALSE)$count,
                     "Proportion" = 0:100,
                     "NbReads" = l)
        }
      }})
    rm(windows)
    null <- sapply(histo,is.null)
    histo <- histo[!null]
    histo <- do.call(rbind,histo)
    if (facet_wrap_chromosomes){
      g <- ggplot2::ggplot(histo, ggplot2::aes(Proportion, Count, fill=NbReads, width=1))+
        scale_fill_discrete(breaks=leg)+
        ggplot2::facet_wrap(~Chr)
    }
    else{
      g <- ggplot2::ggplot(histo, ggplot2::aes(Proportion, Count, fill=NbReads, width=1))+
        scale_fill_discrete(breaks=leg)
    }
    g <- g + ggplot2::geom_bar(stat="identity", width=.3,position="stack") +
      ggplot2::ylab("Count") +
      ggplot2::xlab("Positive Proportion")+
      ggplot2::theme_bw()
    if (save==TRUE){
      message("The plot will be saved to the file ",file)
       ggplot2::ggsave(filename = file,plot = g)
    }
    else{
      g
    }
  }
}
