#include <Rcpp.h>
#include "random"
using namespace Rcpp;

//' @title  Calculate the reads to be kept.
//'
//' @description Calculate the reads to be kept based of the strand proportion of kept windows. This function is called by the functions filterOne or filterMulti.
//' If a kept window has more positive alignments than negative alignments, then 
//' each postive alignment in that window is kept with a probability equal to (positive proportion -  negative proportion)/(positive proportion); and each negative alignment of that window is kept with a probability equal to the given errorRate (0.01 by default). Similarly for windows that have negative alignments more than positive alingments.
//' Since each alignment can be belonged to several windows, then the final probability for keeping of each alignment is the maximum one from all windows containing it.
//' @param posFragments the data frame contains the information of positive fragments (generated from the function getFragment)
//' @param negFragments the data frame contains the information of negative fragments (generated from the function getFragment)
//' @param keptPosWin the data frame contains the information of positive kept windows 
//' @param keptNegWin the data frame contains the information of negative kept windows 
//' @param win the size of the sliding windows
//' @param step the size of the step of sliding windows
//' @param errorRate the probability that an RNA read takes the false strand
//'
//' @return A list of two vectors containing the positive and negative alignments to be kept.
//' 
//' @seealso getFragment, filterOne, filterMulti
//' 
//' @examples
//' posFragments <- data.frame("group"=sample(1:800,1000,replace=TRUE),"start"=sample(1:1000000,1000,replace=TRUE)) 
//' posFragments <- posFragments[order(posFragments$start),] %>% dplyr::mutate(end=start+sample(50:100,1000,replace=TRUE))
//' negFragments <- data.frame("group"=sample(1:800,1000,replace=TRUE),"start"=sample(1:1000000,1000,replace=TRUE)) 
//' negFragments <- negFragments[order(negFragments$start),] %>% dplyr::mutate(end=start+sample(50:100,1000,replace=TRUE))
//' keptPosWin <- data.frame("win"=sample(1:10000,5000,replace=FALSE),"propor"=runif(5000, min=0.7, max=1))
//' keptPosWin <- keptPosWin[order(keptPosWin$win),]
//' keptNegWin <- data.frame("win"=sample(1:10000,5000,replace=FALSE),"propor"=runif(5000, min=0, max=0.3))
//' keptNegWin <- keptNegWin[order(keptNegWin$win),]
//' win <- 1000
//' step <- 100
//' errorRate <- 0.01
//' reads <- rnaCleanR::keepRead(posFragments,negFragments,keptPosWin,keptNegWin,win,step,errorRate)
//' 
//' @export
//' 
// [[Rcpp::export]]


List keepRead(DataFrame posFragments,DataFrame negFragments,DataFrame keptPosWin,DataFrame keptNegWin,int win,int step,double errorRate){
//List keepFragment(IntegerVector startPosAl,IntegerVector endPosAl,IntegerVector groupPos,IntegerVector startNegAl,IntegerVector endNegAl,IntegerVector groupNeg,NumericVector proporWinPos,NumericVector proporWinNeg,IntegerVector keepWinPos,IntegerVector keepWinNeg,int nbWin,int win,int step,double errorRate){
  //Get columns of each data frame
  IntegerVector startPosFragments = posFragments["start"];
  IntegerVector endPosFragments = posFragments["end"];
  IntegerVector groupPosFragments = posFragments["group"];
  IntegerVector startNegFragments = negFragments["start"];
  IntegerVector endNegFragments = negFragments["end"];
  IntegerVector groupNegFragments = negFragments["group"];
  IntegerVector startPos = keptPosWin["Start"];
  NumericVector proporPos = keptPosWin["propor"];
  IntegerVector startNeg = keptNegWin["Start"];
  NumericVector proporNeg = keptNegWin["propor"];
  int nbWin = 0; //number of considering windows
  if (startPos.size()>0) nbWin = startPos[startPos.size()-1];
  if (startNeg.size()>0) nbWin = std::max(nbWin,startNeg[startNeg.size()-1]);
  int* keepWin = new int[nbWin];
  for (int i=0;i<nbWin;i++) keepWin[i]=0;
  for (int i=0;i<keptPosWin.nrows();i++){
    keepWin[startPos[i]-1]=1;
  }
  for (int i=0;i<keptNegWin.nrows();i++){
    if (keepWin[startNeg[i]-1]==1){
      keepWin[startNeg[i]-1]=2;
    }
    else keepWin[startNeg[i]-1]=-1;
  }
  double* probPosWin = new double[nbWin];//probability of positive win
  double* probNegWin = new double[nbWin];//probability of negative win
  for (int i=0;i<nbWin;i++){
    probPosWin[i] = 0;
    probNegWin[i] = 0;
  }
  for (int w=0;w<keptPosWin.nrows();w++){
    int i = startPos[w];
    if (keepWin[i-1]==1){
      probPosWin[i-1]=(2-1/proporPos[w])-errorRate;
      probNegWin[i-1]=errorRate;
    }
    else if (keepWin[i-1]==2){
      probPosWin[i-1]=(1-errorRate);
      probNegWin[i-1]=(1-errorRate);
    }
  }
  for (int w=0;w<keptNegWin.nrows();w++){
    int i = startNeg[w];
    if (keepWin[i-1]==-1){
      probNegWin[i-1]=((1-2*proporNeg[w])/(1-proporNeg[w]))-errorRate;
      probPosWin[i-1]=errorRate;
    }
  }
  double* probPosAl = new double[posFragments.nrows()+1];  
  for (int i=0;i<posFragments.nrows();i++) probPosAl[i]=0;
  for (int a=0;a<posFragments.nrows();a++){
    int s = startPosFragments[a];
    int e = endPosFragments[a];
    int x1 = 1+ceil((double)(s-win)/step);//the first windows that the start position of the read is in
    int x2 = 1+floor((double)(e-1)/step);//the last windows that the end position of the read is in
    if (x1<0) x1=0;
    if (x2>=nbWin) x2=nbWin-1;
    double prob = 0;
    for (int k=x1;k<=x2;k++){ 
      prob=std::max(prob,probPosWin[k]);
    }
    probPosAl[groupPosFragments[a]]=std::max(probPosAl[groupPosFragments[a]],prob);
  }
  double* probNegAl = new double[negFragments.nrows()+1]; 
  for (int i=0;i<negFragments.nrows();i++) probNegAl[i]=0;
  for (int a=0;a<negFragments.nrows();a++){
    int s = startNegFragments[a];
    int e = endNegFragments[a];
    int x1 = 1+ceil((double)(s-win)/step);//the first windows that the start position of the read is in
    int x2 = 1+floor((double)(e-1)/step);//the last windows that the end position of the read is in
    if (x1<0) x1=0;
    if (x2>=nbWin) x2=nbWin-1;
    double prob = 0;
    for (int k=x1;k<=x2;k++){ 
      prob=std::max(prob,probNegWin[k]);
    }
    probNegAl[groupNegFragments[a]]=std::max(probNegAl[groupNegFragments[a]],prob);
  }
  
  std::default_random_engine generator;
  std::vector<int> listKeepPosReads;
  std::vector<int> listKeepNegReads;
  for (int i=0;i<posFragments.nrows();i++){
    std::binomial_distribution<int> distribution(1,probPosAl[groupPosFragments[i]]);
    if (distribution(generator)){
      listKeepPosReads.push_back(groupPosFragments[i]);
    }
  }
  for (int i=0;i<negFragments.nrows();i++){
    std::binomial_distribution<int> distribution(1,probNegAl[groupNegFragments[i]]);
    if (distribution(generator)){
      listKeepNegReads.push_back(groupNegFragments[i]);
    }
  }
  delete[] keepWin;
  delete[] probPosWin;
  delete[] probNegWin;
  delete[] probPosAl;
  delete[] probNegAl;
  return List::create(
    _["Pos"] = listKeepPosReads,
    _["Neg"] = listKeepNegReads);
}
