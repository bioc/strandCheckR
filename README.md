[![Build Status](https://travis-ci.org/UofABioinformaticsHub/strandCheckR.svg?branch=master)](https://travis-ci.org/UofABioinformaticsHub/strandCheckR)

# strandCheckR

This package aims to check the strandeness of the reads in a bam file. It can be applied to quantify and remove putative double strand DNA from a strand-specific RNA sample.
The package uses a sliding window to scan a bam file to get the number of positive/negative reads in each window. 
It then provides method to plot the positive/negative proportions of all sliding windows, which allow users to have an idea about how much the sample was contaminated and the appropriate threshold to be used for filtering. Finally, the users can filter putative DNA contamination from their strand-specific RNA sample using their selected threshold.


# Installation

To install the release version on Bioconductor:
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install('strandCheckR')
```

To install the development version on github:
```
devtools::install_github('UofABioinformaticsHub/strandCheckR', build_vignettes = TRUE)
```

# Usage


