getWin <- function(range,win,step){
  s<-start(range)
  e<-end(range)
  lapply(1:length(s),function(i){
      x1 = 1+ceiling((s[i]-win)/step) #the first windows containing start position 
      x2 = 1+floor((e[i]-1)/step) #the last windows containing the end position
  }) %>% unlist() 
}