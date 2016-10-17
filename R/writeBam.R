writeBam <- function(keep,bamfilein,bamfileout,chromosomes,allChromosomes,lenSeq){#write the filter bamfile based on the selected keep reads
  for (i in c(1:length(bamfilein))){
    reader <- bamReader(bamfilein[i],idx=TRUE) #open a reader of the input bamfile to extract read afterward
    header <- getHeader(reader) #get the header of the input bam file
    writer <- bamWriter(header,bamfileout[i]) #prepare to write the output bamfile with the same header
    message("Writing sample ",i)
    for (chr in chromosomes){
      chromosomeIndex <- which(allChromosomes==chr)
      if (length(keep[[chr]][[i]])>0){
        #get the range of kept reads
        range <- bamRange(reader,c(chromosomeIndex-1,0,lenSeq[chromosomeIndex]))
        #write the kept reads into output file
        bamSave(writer,range[keep[[chr]][[i]],],refid=chromosomeIndex-1)
        keep[[chr]][[i]] <- c(1)
        remove(range)
      }
    }
    bamClose(writer)
    remove(reader)
    remove(header)
  }
  remove(keep)
}
