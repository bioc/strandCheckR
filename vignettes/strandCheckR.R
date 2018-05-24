## ----getWin, eval=FALSE----------------------------------------------------
#  library(strandCheckR)
#  file <- system.file("extdata","s1.sorted.bam",package = "strandCheckR")
#  win <- getWinFromBamFile(file)

## ----plot, eval=FALSE------------------------------------------------------
#  plotHist(win,save=TRUE,file="hist.png",width=8,height=5)
#  plotWin(win,save=TRUE,file="win.png",width=8,height=5)

## ----filteDNA, eval=FALSE--------------------------------------------------
#  filterDNA(bamfilein = file,bamfileout = "filter.bam", threshold = 0.7)

