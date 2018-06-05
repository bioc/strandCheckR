## ----getWin, message=FALSE,warning=FALSE-----------------------------------
library(strandCheckR)
files <- system.file("extdata",c("s1.sorted.bam","s2.sorted.bam"),
                    package = "strandCheckR")
win <- getWinFromBamFile(files)
# shorten the file name
win$File <- basename(win$File)
win

## ----intersect, eval=TRUE, warning=FALSE,message=FALSE---------------------
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
annot <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
annot <- annot[!duplicated(annot),]
#add chr before the chromosome names to be consistent with the annot data
win$Seq <- paste0("chr",win$Seq) 
win <- intersectWithFeature(windows = win, annotation = annot, 
                            overlapCol = "OverlapTranscript")
win

## ----intersectGetFeature, eval=TRUE, warning=FALSE,message=FALSE-----------
win <- intersectWithFeature(windows = win, annotation = annot,
                            mcolsAnnot = "tx_name",collapse = ",")
win[order(win$MaxCoverage,decreasing = TRUE),
    c("Seq","Start","End","MaxCoverage","File","tx_name")]

## ----plotHist, eval=TRUE, message=FALSE,warning=FALSE----------------------
hist <- summarizeHist(windows = win,group_by = c("File","OverlapTranscript"), 
                    normalize_by = "File")
plotHist(hist, facets = c("File","OverlapTranscript"), scales = "free_y")

## ----heatMap, eval=TRUE,fig.height = 6, fig.width = 14,fig.align = "center"----
plotHist(hist, facets = c("OverlapTranscript"), heatmap = TRUE)

## ----plotwin,eval=TRUE,message=FALSE,warning=FALSE-------------------------
plotWin(win, facets = c("File","OverlapTranscript"))

## ----filterDNA, eval=TRUE, message=FALSE, warning=FALSE, results=FALSE-----
filterDNA(file = files[1], destination = "s1.filter.bam", threshold = 0.7)

