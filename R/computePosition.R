computePosition <- function(alignments,chr){
  positionPos <- data.frame("group"=c(),"start"=c(),"end"=c(),"sample"=c()) #position of each fragment of positive reads
  positionNeg <- data.frame("group"=c(),"start"=c(),"end"=c(),"sample"=c()) #position of each fragment of negative reads
  index <- list() #the index of each read in the alignment from its index in positive/negative read lists
  nbRead <- c() #the number of reads
  for (i in c(1:length(alignments))){ #compute positionPos and positionNeg for all samples
    al <- alignments[[i]][seqnames(alignments[[i]])==chr]
    nbRead <- c(nbRead,length(al))
    index[[i]] <- getIndex(as.vector(strand(al)))
    alPos<-al[strand(al)=="+"]
    positionPos <- rbind(positionPos,extractAlignmentRangesOnReference(cigar(alPos),pos=start(alPos)) 
                         %>% data.frame() 
                         %>% dplyr::select(-c(group_name,width)) 
                         %>% mutate(sample=i))
    remove(alPos)
    gc()
    alNeg<-al[strand(al)=="-"]
    remove(al)
    gc()
    positionNeg <- rbind(positionNeg,extractAlignmentRangesOnReference(cigar(alNeg),pos=start(alNeg)) 
                         %>% data.frame() 
                         %>% dplyr::select(-c(group_name,width))
                         %>% mutate(sample=i))
    remove(alNeg)
    gc()
  } 
  return (list(nbRead=nbRead,Pos=positionPos,Neg=positionNeg,Index=index))
}
