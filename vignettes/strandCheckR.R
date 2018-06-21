## ----getWin, message=FALSE,warning=FALSE-----------------------------------
library(strandCheckR)
files <- system.file("extdata",c("s1.sorted.bam","s2.sorted.bam"),
                    package = "strandCheckR")
win <- getWinFromBamFile(files, sequences = "10")
# shorten the file name
win$File <- basename(win$File)
win

## ----highestCoverage, eval=TRUE, message=FALSE,warning=FALSE---------------
head(win[order((win$NbPositive+win$NbNegative),decreasing=TRUE),])

## ----paired end, eval=TRUE,message=FALSE,warning=FALSE---------------------
fileP <- system.file("extdata","paired.bam",package = "strandCheckR")
winP <- getWinFromBamFile(fileP, sequences = "10")
winP$File <- basename(winP$File) #shorten file name
winP

## ----intersect, eval=TRUE, warning=FALSE,message=FALSE---------------------
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
annot <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
#add chr before the sequence names to be consistent with the annot data
win$Seq <- paste0("chr",win$Seq) 
win <- intersectWithFeature(windows = win, annotation = annot, 
                            overlapCol = "OverlapTranscript")
win

## ----plotHist, eval=TRUE, message=FALSE,warning=FALSE----------------------
plotHist(windows = win, group_by = c("File","OverlapTranscript"), 
        normalize_by = "File", scales = "free_y")

## ----plotHistPaired, eval=TRUE,message=FALSE,warning=FALSE-----------------
plotHist(windows = winP, group_by = "Type", normalize_by = "Type", 
        scales = "free_y")

## ----heatMap, eval=TRUE, message = FALSE, warning=FALSE--------------------
plotHist(windows = win, group_by = c("File","OverlapTranscript"), 
        normalize_by = "File", heatmap = TRUE)

## ----plotwin,eval=TRUE,message=FALSE,warning=FALSE-------------------------
plotWin(win, group_by = "File")

## ----filterDNA, eval=TRUE, message=FALSE, warning=FALSE, results=FALSE-----
win2 <- filterDNA(file = files[2], sequences = "10", 
        destination = "s2.filter.bam", threshold = 0.7, getWin = TRUE)

## ----compare,eval=TRUE,message=FALSE,warning=FALSE-------------------------
win2$File <- basename(win2$File)
win2$File <- factor(win2$File,levels = c("s2.sorted.bam","s2.filter.bam"))
library(ggplot2)
plotHist(win2,group_by = "File",normalize_by = "File",scales = "free_y") +
    ggtitle("Histogram of positive proportions over all sliding windows before 
            and after filtering reads coming from double strand DNA")

