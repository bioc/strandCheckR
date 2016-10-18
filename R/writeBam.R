writeBam <- function(keep,bamfilein,bamfileout,chromosomes,allChromosomes,lenSeq) {
  #write the filter bamfile based on the selected keep reads
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
      keep[[chr]] <- c(1)
      remove(range)
    }
  }
  bamClose(writer)
  bamClose(reader)
  remove(header)
}
