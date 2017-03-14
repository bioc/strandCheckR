## ----getplot,  eval=FALSE------------------------------------------------
#  file <- system.file("data","s1.chr1.bam",package = "rnaCleanR")
#  getPlot(bamfilein = file,readLength = 100,histPlotFile = "hist.pdf",winPlotFile = "win.pdf")

## ----filterOne, eval=FALSE-----------------------------------------------
#  file <- system.file("data","s1.chr1.bam",package = "rnaCleanR")
#  filterOne(bamfilein = file,bamfileout = "filter.bam", readLength = 100, threshold = 0.7)

## ----filterOnePlot, eval=FALSE-------------------------------------------
#  filterOne(bamfilein = file,bamfileout = "filter.bam", readLength = 100, threshold = 0.7,histPlot = TRUE,histPlotFile = "hist.pdf",winPlot = TRUE,winPlotFile = "win.pdf")

## ----filterMulti, eval=FALSE---------------------------------------------
#  file <- system.file("data",c("s1.chr1.bam","s2.chr1.bam"),package = "rnaCleanR")
#  filterMulti(bamfilein = file,bamfileout = c("filter1.bam","filter2.bam"), readLength = 100, threshold = 0.7)
#  

