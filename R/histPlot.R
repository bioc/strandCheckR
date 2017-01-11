histPlot <- function(windows,histfile,breaks){
  allCols <- c("blue","green","cyan","pink","orange","red","brown","black")
  leg <- c("0-10","10-20","20-50","50-100","100-200","200-500","500-1000",">1000")
  windows$propor <- round(windows$propor*breaks)
  histo <- data.frame("propor"=c(),"group"=c(),"freq"=c())
  for (i in 0:breaks){
    a <- dplyr::filter(windows,propor==i) 
    if (nrow(a)>0){
      a <- table(a$group) %>% as.data.frame()
      histo <- rbind(histo,data.frame(propor = i,group=a$Var1,freq=a$Freq))  
    }
  }
  presentGroup <-  as.numeric(levels(unique(histo$group)))
  allCols <- allCols[presentGroup]
  leg <- leg[presentGroup]
  
  ggplot2::ggplot(histo, ggplot2::aes(propor, freq, fill=group, width=0.9))+ggplot2::geom_bar(stat="identity", width=.3,position="stack") + ggplot2::scale_fill_manual(values =  allCols,labels = leg) + ggplot2::labs(fill = "max Coverage") + ggplot2::ylab("Positive Proportion") + ggplot2::xlab("Frequency")
  ggplot2::ggsave(histfile)#,width=1500,height=700)
}