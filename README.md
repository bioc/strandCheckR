
<!-- README.md is generated from README.Rmd. Please edit the Rmd file only -->
[![Build Status](https://travis-ci.org/UofABioinformaticsHub/strandCheckR.svg?branch=master)](https://travis-ci.org/UofABioinformaticsHub/strandCheckR) [![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![DOI](https://zenodo.org/badge/70646093.svg)](https://zenodo.org/badge/latestdoi/70646093) [![DOI](http://joss.theoj.org/papers/10.21105/joss.01145/status.svg)](https://doi.org/10.21105/joss.01145)

strandCheckR
------------

This package aims to check the strandedness of reads in a bam file, enabling easy detection of any contaminating genomic DNA or other unexpected sources of contamination. It can be applied to quantify and remove reads which correspond to putative double strand DNA within a strand-specific RNA sample. The package uses a sliding window to scan a bam file and find the number of positive/negative reads in each window. It then provides method to plot the proportions of positive/negative stranded alignments within all windows, which allow users to determine how much the sample was contaminated, and to determine an appropriate threshold for filtering. Finally, users can filter putative DNA contamination from any strand-specific RNAseq sample using their selected threshold.

Installation
------------

To install the release version from Bioconductor:

``` r
install.packages("BiocManager")
BiocManager::install("strandCheckR")
```

To install the development version on github (i.e. this version):

``` r
install.packages("BiocManager")
BiocManager::install("UofABioinformaticsHub/strandCheckR")
```

Quick Usage Guide
-----------------

Following are the main functions of the package.

-   `getStrandFromBamFile()`

To get the number of +/- stranded reads of all sliding windows across a bam file:

``` r
# Load the package and example bam files
library(strandCheckR)
files <- system.file(
    "extdata", c("s1.sorted.bam", "s2.sorted.bam"),
    package = "strandCheckR"
)

# Find the read proportions from chromosome 10 for the two files
win <- getStrandFromBamFile(files, sequences = "10")

# Tidy up the file name for prettier output
win$File <- basename(as.character(win$File))
win
## DataFrame with 3078 rows and 10 columns
##       Type   Seq     Start       End NbPos NbNeg CovPos CovNeg MaxCoverage
##      <Rle> <Rle> <numeric> <numeric> <Rle> <Rle>  <Rle>  <Rle>       <Rle>
## 1       SE    10   7696701   7697700     0    17      0    393          17
## 2       SE    10   7696801   7697800     0    17      0    393          17
## 3       SE    10   7696901   7697900     0    17      0    393          17
## 4       SE    10   7697001   7698000     0    17      0    393          17
## 5       SE    10   7697101   7698100     0    17      0    393          17
## ...    ...   ...       ...       ...   ...   ...    ...    ...         ...
## 3074    SE    10   7398501   7399500    46    34   2241   1668          13
## 3075    SE    10   7398601   7399600    46    34   2241   1668          13
## 3076    SE    10   7398701   7399700    41    32   2046   1568          13
## 3077    SE    10   7398801   7399800    48    31   2500   1681          25
## 3078    SE    10   7398901   7399900    52    35   2581   1728          25
##               File
##        <character>
## 1    s1.sorted.bam
## 2    s1.sorted.bam
## 3    s1.sorted.bam
## 4    s1.sorted.bam
## 5    s1.sorted.bam
## ...            ...
## 3074 s2.sorted.bam
## 3075 s2.sorted.bam
## 3076 s2.sorted.bam
## 3077 s2.sorted.bam
## 3078 s2.sorted.bam
```

-   `plotHist()`

The histogram plot shows you the proportion of +/- stranded reads across all windows.

``` r
plotHist(
        windows = win, 
        groupBy = "File", 
        normalizeBy = "File", 
        scales = "free_y"
        )
```

![](README_files/figure-markdown_github/plotHist-1.png)

In this example, *s2.sorted.bam* seems to be contaminated with double stranded DNA, as evidenced by many windows containing a roughly equal proportion of reads on both strands, whilst *s1.sorted.bam* is cleaner.

-   `plotWin()`

The output from `plotWin()` represents each window as a point. This plot also has threshold lines which can be used to provide guidance as to the best threshold to choose when filtering windows. Given a suitable threshold, reads from a positive (resp. negative) window are kept if and only if the proportion is above (resp. below) the corresponding threshold line.

``` r
plotWin(win, groupBy = "File")
```

![](README_files/figure-markdown_github/plotWin-1.png)

-   `filterDNA()`

The function `filterDNA()` removes potential double stranded DNA from a bam file using a selected threshold.

``` r
win2 <- filterDNA(
    file = files[2], 
    destination = "s2.filtered.bam", 
    sequences = "10", 
    threshold = 0.7, 
    getWin = TRUE
)
```

Comparing the histogram plot of the file before and after filtering shows that reads from the windows with roughly equal proportions of +/- stranded reads have been removed.

``` r
win2$File <- basename(as.character(win2$File))
win2$File <- factor(win2$File, levels = c("s2.sorted.bam", "s2.filtered.bam"))
library(ggplot2)
plotHist(win2, groupBy = "File", normalizeBy = "File", scales = "free_y") 
```

![](README_files/figure-markdown_github/plotHistAfterFilter-1.png)

A more comprehensive vignette is available at <https://bioconductor.org/packages/release/bioc/vignettes/strandCheckR/inst/doc/strandCheckR.html>

Support
-------

We recommend that questions seeking support in using the software are posted to the Bioconductor support forum - <https://support.bioconductor.org/> - where they will attract not only our attention but that of the wider Bioconductor community.

Code contributions, bug reports and feature requests are most welcome. Please make any pull requests against the master branch at <https://github.com/UofABioinformaticsHub/strandCheckR> and file issues at <https://github.com/UofABioinformaticsHub/strandCheckR/issues>

Author Contributions
--------------------

-   *Thu-Hien To* authored the vast majority of code within the package along with unit tests
-   *Thu-Hien To* and *Stephen Pederson* worked closely together on the package design and methodology

License
-------

`strandCheckR` is licensed under [GPL &gt;= 2.0](https://www.r-project.org/Licenses/GPL-2)
