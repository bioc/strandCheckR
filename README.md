[![Build Status](https://travis-ci.org/UofABioinformaticsHub/strandCheckR.svg?branch=master)](https://travis-ci.org/UofABioinformaticsHub/strandCheckR)

# strandCheckR

An R Package that uses a sliding windows approach to scan a bam file and get the number of positive/negative reads in each window. 

Applicable to check the strandedness of a stranded RNA-Seq experiment and to filter reads coming from potential double strand DNA if there's any. 

# Installation

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(c("devtools", "BiocStyle", "GenomeInfoDb", "GenomicAlignments", "GenomicRanges", "IRanges", "Rsamtools", "S4Vectors", "BiocGenerics", "TxDb.Hsapiens.UCSC.hg38.knownGene"))
devtools::install_github('UofABioinformaticsHub/strandCheckR', build_vignettes = TRUE)
library(strandCheckR)
```



