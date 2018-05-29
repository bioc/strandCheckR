[![Build Status](https://travis-ci.org/UofABioinformaticsHub/strandCheckR.svg?branch=master)](https://travis-ci.org/UofABioinformaticsHub/strandCheckR)

# strandCheckR

An R Package that uses a sliding windows approach to scan a bam file and get the number of positive/negative reads in each window. Applicable to check the strandedness of a stranded RNA-Seq experiment and to filter reads coming from potential double strand DNA if there's any. 

# Installation

```
devtools::install_github('UofABioinformaticsHub/strandCheckR', build_vignettes = TRUE)
library(strandCheckR)
```

# Vignette

The vignette for usage is [here](https://uofabioinformaticshub.github.io/strandCheckR/vignettes/strandCheckR.html)


