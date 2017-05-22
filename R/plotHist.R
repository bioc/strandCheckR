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
  
  if ("Type" %in% colnames(windows)){
    histoFirst <- lapply(leg,function(l){
      a <- dplyr::filter(windows,MaxCoverage==l,Type=="First")$Proportion
      if (length(a)>0){
        hist(a,breaks = 100)$count  
      }}) 
    null <- sapply(histoFirst,is.null)
    histoFirst <- histoFirst[!null] %>% 
      as.data.frame() %>%  
      set_colnames(leg[!null]) %>%
      dplyr::mutate("Type"="First","Proportion"=1:100) %>% 
      melt(id.vars = c("Type","Proportion")) %>% 
      set_colnames(c("Type","Proportion","MaxCoverage","Count"))
    
    histoSecond <- lapply(leg,function(l){
      a <- dplyr::filter(windows,MaxCoverage==l,Type=="Second")$Proportion
      if (length(a)>0){
        hist(a,breaks = 100)$count  
      }}) 
    null <- sapply(histoSecond,is.null)
    histoSecond <- histoSecond[!null] %>% 
      as.data.frame() %>%  
      set_colnames(leg[!null]) %>%
      dplyr::mutate("Type"="Second","Proportion"=1:100) %>% 
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
        hist(a,breaks = 100)$count  
      }}) 
    null <- sapply(histo,is.null)
    histo <- histo[!null] %>% 
      as.data.frame() %>%  
      set_colnames(leg[!null]) %>%
      dplyr::mutate("Proportion"=1:100) %>% 
      reshape2::melt(id.vars = c("Proportion")) %>% 
      set_colnames(c("Proportion","MaxCoverage","Count"))
    ggplot(histo, aes(Proportion, Count, fill=MaxCoverage, width=1))+
      geom_bar(stat="identity", width=.3,position="stack") + 
      ylab("Count") + 
      xlab("Positive Proportion")+ 
      theme_bw() 
  }
}
