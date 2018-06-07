# #' @title get window data frame with the correct sequence name and position
# #' 
# #' @description get the correct sequence name and position for each window
# #' 
# #' @param Win a data frame contains the information of every window in 
# #' \code{sequences}
# #' @param sequences a vector of sequence names
# #' @param sequenceInfo a data frame that contains some information of 
# #' the alignments
# #' @param winWidth the length of sliding window
# #' @param winStep the step length to sliding the window
# #' 
getWinInSequence <- function(Win, sequences, sequenceInfo, winWidth = 1000, 
                            winStep = 100){

    # Check the correct columns are in the sequenceInfo df
    reqCols <- c("FirstBaseInPartition", "LastBaseInPartition")
    if (!all(reqCols %in% names(sequenceInfo))){ 
        stop("sequenceInfo must contain the columns ", reqCols)
    }
    stopifnot(is.numeric(winWidth) || is.numeric(winStep))

    for (i in seq_along(sequences)){
        if (!is.na(sequenceInfo$FirstBaseInPartition[i])){
            currentChr <- sequences[i]
            id <- which(sequences == currentChr)
            # get id of the first window of the sequence
            idFirst <- ceiling(sequenceInfo$FirstBaseInPartition[id] / winStep) 
            # get id of the last window of the sequence
            idLast <- ceiling(
                (sequenceInfo$LastBaseInPartition[id] - winWidth+1) / winStep)
            #get the windows of the sequence
            idRows <- which(Win$Start >= idFirst & Win$Start <= idLast) 
            Win$Seq[idRows] <- Rle(currentChr)
            Win$Start[idRows] <- (Win$Start[idRows] - idFirst)+1
        }
    }
    return(Win[Win$Seq!="",])
}
