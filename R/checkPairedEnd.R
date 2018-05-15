#' @title Test wheather a bam file if single end or paired end
#' @description Check the first 100000 first reads of the bam file to see wheather it is single end or paired end 
#' @param file the input bam file. Your bamfile should be sorted and have an index file located at the same path as well.
#' @param yieldSize the number of reads to be checked, 100000 by default.
#' @export
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools scanBam
#' @importFrom Rsamtools ScanBamParam

checkPairedEnd <- function(file, yieldSize = 100000){
  message("Testing paired end bam file by checking the first ", yieldSize," reads")
  checkFile <- BamFile(file, yieldSize = yieldSize)
  flag <- scanBam(checkFile, param = ScanBamParam(what = "flag"))[[1]]$flag
  paired <- any(flag %% 2 == 1)
  if (paired) message("Your bam file is paired end") else message("Your bam file is single end")  
  return(paired)
}

