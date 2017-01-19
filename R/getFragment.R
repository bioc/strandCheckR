#' @title Generating spliced fragments of the input alignments
#' 
#' @description This method use the cigar and the start mapping position of each reads to generate the position of each spliced alignment

#' @param alignment the alignments to get spliced fragments
#' 
#' @return a list of two data frame for fragments come from positive and negative reads. Each data frame has information
#' about the start/end position of the fragments and the group which allows the refer which read that the fragment comes from
#' 
#' @examples 
#' bamfilein <- system.file("data","s1.chr1.bam",package = "rnaCleanR")
#' alignment <- GenomicAlignments::readGAlignments(bamfilein) 
#' alignmentInChr1 <- alignment[seqnames(alignment)=="1"] 
#' fragments <- getFragment(alignmentInChr1)
#' 
#' @export
#'
getFragment <- function(alignment){
  position <- data.frame("group"=c(),"start"=c(),"end"=c()) #data frame contains position of each fragment of reads
  alPos <- alignment[strand(alignment)=='+']
  alNeg <- alignment[strand(alignment)=='-']
  positionPos <- GenomicAlignments::extractAlignmentRangesOnReference(GenomicAlignments::cigar(alPos),pos=start(alPos)) %>% data.frame() %>% dplyr::select(-c(group_name,width))  
  positionPos <- positionPos[order(positionPos$start),]
  positionNeg <- GenomicAlignments::extractAlignmentRangesOnReference(GenomicAlignments::cigar(alNeg),pos=start(alNeg)) %>% data.frame() %>% dplyr::select(-c(group_name,width))  
  positionNeg <- positionNeg[order(positionNeg$start),]
  return(list(Pos = positionPos,Neg = positionNeg))
}