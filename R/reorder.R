reorder <- function(position,win,step,lim){
  position$Pos[["firstW"]] <- ceiling((position$Pos$start-win-1+floor((position$Pos$end-position$Pos$start+1)*lim))/step)
  position$Neg[["firstW"]] <- ceiling((position$Neg$start-win-1+floor((position$Neg$end-position$Neg$start+1)*lim))/step)
  position$Pos <-  position$Pos[order(position$Pos$firstW),] #reorder positionPos following the starting position of each fragment
  position$Neg <-  position$Neg[order(position$Neg$firstW),] #reorder positionNeg following the starting position of each fragment
  return (position)
}
