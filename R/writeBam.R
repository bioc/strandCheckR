#' @title Write Bam File
#' 
#' @description write the reads into the output bamfile
#' 
#' @details None
#' 
#' @param keep the indices of the reads in the input bamfile to be kept
#' @param bamfilein the input bam file 
#' @param bamfileout the output bam file
#' @param chromosomes the list of chromosomes to write
#' @param allChromosomes the list of all chromosomes
#' @param lenSeq the length of every chromosome
#'
writeBam <- function(keep,bamfilein,bamfileout,chromosomes,allChromosomes,lenSeq) {
  reader <- bamReader(bamfilein, idx = TRUE) #open a reader of the input bamfile to extract read afterward
  header <- getHeader(reader) #get the header of the input bam file
  writer <- bamWriter(header, bamfileout) #prepare to write the output bamfile with the same header
  for (chr in chromosomes) {
    chromosomeIndex <- which(allChromosomes == chr)
    if (length(keep[[chr]]) > 0) {
      #get the range of kept reads
      range <- bamRange(reader, c(chromosomeIndex - 1, 0, lenSeq[chromosomeIndex]))
      #write the kept reads into output file
      bamSave(writer, range[keep[[chr]], ], refid = chromosomeIndex - 1)
      remove(range)
      keep[[chr]] <- c(1)
    }
  }
  bamClose(writer)
  bamClose(reader)
  rm(list=ls())
  gc()
}
