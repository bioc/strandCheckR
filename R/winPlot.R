winPlot <- function(windows,winfile,xlim){
  allCols <- c("blue","green","cyan","pink","orange","red","brown","black")
  leg <- c("0-10","10-20","20-50","50-100","100-200","200-500","500-1000",">1000")
  windows$propor <- round(windows$propor,2)
  windows$sum <- round(windows$sum)
  windows <- windows[!duplicated(windows),]
  pdf(winfile)
  m <- max(windows$sum)
  if (!missing(xlim) && xlim < m){
    m <- xlim
    windows <- dplyr::filter(windows,sum<=m)
  }
  plot(x=windows$sum,y = windows$propor,ylim = c(0,1),xlim = c(1,m),col = allCols[windows$group],pch=16,xlab="Sum",ylab="Positive Proportion")
  n <- seq(m)
  colors <- c("#00441b","#238b45","#74c476","#c7e9c0")
  i <- 1
  for (p in c(0.6,0.7,0.8,0.9)){
    se <- sqrt(p*(1-p)/n)
    upper <- p + qnorm(0.975)*se
    lines(n, upper, type = "l", ylim = c(0, 1),xlim = c(1,m), col = colors[i])
    lines(n, 1-upper, type = "l", ylim = c(0, 1),xlim = c(1,m), col = colors[i])
    text(m-20,p+0.01,labels = as.character(p))
    text(m-20,1-p,labels = as.character(p))
    i <- i+1
  }
  legend("center",title = "max Coverage", col=allCols,leg=leg,pt.cex=1,cex=1,pch = 16)
  dev.off()
}

