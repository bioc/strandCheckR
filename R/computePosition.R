computePosition <- function(alignments){
  positionPos <- data.frame("group"=c(),"start"=c(),"end"=c(),"sample"=c()) #position of each fragment of positive reads
  positionNeg <- data.frame("group"=c(),"start"=c(),"end"=c(),"sample"=c()) #position of each fragment of negative reads
  index <- list() #the index of each read in the alignment from its index in positive/negative read lists
  nbRead <- c() #the number of reads
  for (i in c(1:length(alignments))){ #compute positionPos and positionNeg for all samples
    nbRead <- c(nbRead,length(alignments[[i]]))
    index[[i]] <- getIndex(as.vector(strand(alignments[[i]])))
    alPos<-alignments[[i]][strand(alignments[[i]])=="+"]
    positionPos <- rbind(positionPos,extractAlignmentRangesOnReference(cigar(alPos),pos=start(alPos)) 
                         %>% data.frame() 
                         %>% dplyr::select(-c(group_name,width)) 
                         %>% mutate(sample=i))
    remove(alPos)
    alNeg<-alignments[[i]][strand(alignments[[i]])=="-"]
    positionNeg <- rbind(positionNeg,extractAlignmentRangesOnReference(cigar(alNeg),pos=start(alNeg)) 
                         %>% data.frame() 
                         %>% dplyr::select(-c(group_name,width))
                         %>% mutate(sample=i))
    remove(alNeg)
  } 
  rm(list=setdiff(ls(), c("nbRead","positionPos","positionNeg","index")))
  gc()
  return (list(nbRead=nbRead,Pos=positionPos,Neg=positionNeg,Index=index))
}
