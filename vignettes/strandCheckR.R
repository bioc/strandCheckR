## ----getWin, message=FALSE,warning=FALSE-----------------------------------
library(strandCheckR)
files <- system.file("extdata",c("s1.sorted.bam","s2.sorted.bam"),
                    package = "strandCheckR")
win <- getWinFromBamFile(files)
# shorten the file name
win$File <- basename(win$File)
win

## ----highestCoverage, eval=TRUE, message=FALSE,warning=FALSE---------------
head(win[order((win$NbPositive+win$NbNegative),decreasing=TRUE),])

## ----intersect, eval=TRUE, warning=FALSE,message=FALSE---------------------
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
annot <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
annot <- annot[!duplicated(annot),]
#add chr before the sequence names to be consistent with the annot data
win$Seq <- paste0("chr",win$Seq) 
win <- intersectWithFeature(windows = win, annotation = annot, 
                            overlapCol = "OverlapTranscript")
win

## ----intersectGetFeature, eval=TRUE, warning=FALSE,message=FALSE-----------
win <- intersectWithFeature(windows = win, annotation = annot,
                            mcolsAnnot = "tx_name",collapse = ",")
winOverlap <- win[win$tx_name!="unknown",]
winOverlap[order(winOverlap$MaxCoverage,decreasing = TRUE),
    c("Seq","Start","End","MaxCoverage","File","tx_name")] 

## ----plotHist, eval=TRUE, message=FALSE,warning=FALSE----------------------
hist <- summarizeHist(windows = win,group_by = c("File","OverlapTranscript"), 
                    normalize_by = "File")
plotHist(hist, facets = c("File","OverlapTranscript"), scales = "free_y")

## ----heatMap, eval=TRUE, fig.width = 25, fig.height=10, warning=FALSE------
plotHist(hist, facets = c("OverlapTranscript"), heatmap = TRUE)

## ----plotwin,eval=TRUE,message=FALSE,warning=FALSE-------------------------
plotWin(win, facets = "File")

## ----filterDNA, eval=TRUE, message=FALSE, warning=FALSE, results=FALSE-----
filterDNA(file = files[1], destination = "s1.filter.bam", threshold = 0.7)

