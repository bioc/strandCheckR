getStart <- function(range,win,step){
  s<-start(range)
  e<-end(range)
  lapply(1:length(s),function(i){
        x1 = 1+step*ceiling((s[i]-win)/step) #the starting base of the first windows containing s
        x2 = 1+step*floor((e[i]-1)/step) #the starting base of the last windows containing s
      }) %>% unlist() 
} 