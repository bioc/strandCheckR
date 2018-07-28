# #' @title get window data frame with the correct sequence name and position 
# #' @description get the correct sequence name and position for each window 
# #' @param Win a data frame contains the information of every window in 
# #' \code{sequences} 
# #' @param seqInfo a data frame that contains some information of 
# #' the alignments
# #' @param winWidth the length of sliding window 
# #' @param winStep the step
# #' length to sliding the window 
getWinInSequence <- function(Win, seqInfo, winWidth = 1000, winStep = 100) 
{
    # Check the correct columns are in the seqInfo df
    reqCols <- c("FirstBaseInPart", "LastBaseInPart")
    if (!all(reqCols %in% names(seqInfo))) {
        stop("seqInfo must contain the columns ", reqCols)
    }
    stopifnot(is.numeric(winWidth) || is.numeric(winStep))
    
    for (id in seq_along(seqInfo$Sequence)) {
        currentChr <- seqInfo$Sequence[id]
        # get id of the first window of the sequence
        idFirst <- ceiling(seqInfo$FirstBaseInPart[id]/winStep)
        # get id of the last window of the sequence
        idLast <- ceiling((seqInfo$LastBaseInPart[id] - winWidth + 1)/winStep)
        # get the windows of the sequence
        idRows <- which(Win$Start >= idFirst & Win$Start <= idLast)
        Win$Seq[idRows] <- Rle(currentChr)
        Win$Start[idRows] <- (Win$Start[idRows] - idFirst) + 1
    }
    return(Win[Win$Seq!="",])
}
