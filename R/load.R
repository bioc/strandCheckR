.onLoad <- function(libname=find.package("strandCheckR"), pkgname = "strandCheckR"){
  utils::globalVariables(c("group_name",
                           "CovPositive","CovNegative","NbPositive","NbNegative","MaxCoverage","Type","PositiveProportion","Count",
                           "Chr","NbReads","Threshold","ThresholdP","ThresholdN","PositiveProportion")) 
  invisible()
}