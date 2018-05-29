## ----getWin, eval=FALSE----------------------------------------------------
#  library(strandCheckR)
#  files <- system.file("extdata",c("s1.sorted.bam","s2.sorted.bam"),package = "strandCheckR")
#  win <- getWinFromBamFile(files)
#  win$File <- basename(win$File)

## ----intersect, eval=FALSE-------------------------------------------------
#  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#  annot <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
#  win$Chr <- paste0("chr",win$Chr) #add chr before the chromosome names to be consistent with the annot data
#  win <- intersectWithFeature(windows = win,annotation = annot,overlapCol = "OverlapFeature")

## ----plotHist, eval=FALSE--------------------------------------------------
#  hist <- summarizeHist(windows = win,group_by = c("File","OverlapFeature"), normalize_by = "File")
#  plotHist(hist, facets = c("File","OverlapFeature"), save=TRUE,file="hist.png",width=10,height=7)

## ----heatMap,eval=FALSE----------------------------------------------------
#  plotHist(hist, facets = c("OverlapFeature"), heatmap = TRUE, save=TRUE,file="histHeat.png",width=16,height=8)

## ----plotwin,eval=FALSE----------------------------------------------------
#  plotWin(win, facets = c("File","OverlapFeature"), save=TRUE,file="win.png",width=10,height=7)

## ----filterDNA, eval=FALSE-------------------------------------------------
#  filterDNA(file = files[1], destination = "s1.filter.bam", threshold = 0.7)

