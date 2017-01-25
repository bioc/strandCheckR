#include <Rcpp.h>
using namespace Rcpp;

//' @title  Calculate the reads to be kept based of the kept windows
//'
//' @description 
//'
//' @param posFragments the data frame contains the information of positive fragments (generated from the function getFragment)
//' @param negFragments the data frame contains the information of negative fragments (generated from the function getFragment)
//' @param keptPosWin the data frame contains the information of windows that we'll keep positive fragments
//' @param keptNegWin the data frame contains the information of windows that we'll keep negative fragments
//' @param win the size of the sliding windows
//' @param step the size of the step of sliding windows
//' @param errorRate the probability that an RNA read takes the false strand
//'
//' @return A list of two vectors containing the positive and negative reads to be kept
//' 
//' @examples
//' posFragments <- data.frame("group"=sample(1:800,1000,replace=TRUE),"start"=sample(1:1000000,1000,replace=TRUE)) %>% dplyr::mutate(end=start+sample(50:100,1000,replace=TRUE))
//' negFragments <- data.frame("group"=sample(1:800,1000,replace=TRUE),"start"=sample(1:1000000,1000,replace=TRUE)) %>% dplyr::mutate(end=start+sample(50:100,1000,replace=TRUE))
//' keptPosWin <- data.frame("win"=sample(1:10000,5000,replace=TRUE),"propor"=runif(5000, min=0, max=1))
//' keptNegWin <- data.frame("win"=sample(1:10000,5000,replace=TRUE),"propor"=runif(5000, min=0, max=1))
//' win <- 1000
//' step <- 100
//' errorRate <- 0.01
//' keepRead(posFragments,negFragments,keptPosWin,keptNegWin,win,step,errorRate)
//' 
//' @export
// [[Rcpp::plugins(cpp11)]]
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
  IntegerVector winPos = keptPosWin["win"];
  NumericVector proporPos = keptPosWin["propor"];
  IntegerVector winNeg = keptNegWin["win"];
  NumericVector proporNeg = keptNegWin["propor"];
  int nbWin = 0; //number of considering windows
  if (winPos.size()>0) nbWin = winPos[winPos.size()-1];
  if (winNeg.size()>0) nbWin = std::max(nbWin,winNeg[winNeg.size()-1]);
    
  int* keepWin = new int[nbWin];
  for (int i=0;i<nbWin;i++) keepWin[i]=0;
  for (int i=0;i<keptPosWin.nrows();i++){
    keepWin[winPos[i]-1]=1;
  }
  for (int i=0;i<keptNegWin.nrows();i++){
    if (keepWin[winNeg[i]-1]==1){
      keepWin[winNeg[i]-1]=2;
    }
    else keepWin[winNeg[i]-1]=-1;
  }
  NumericVector probPosWin(nbWin);//probability of positive win
  NumericVector probNegWin(nbWin);//probability of negative win
  for (int i=0;i<nbWin;i++){
    probPosWin[i] = 0;
    probNegWin[i] = 0;
  }
  for (int w=0;w<keptPosWin.nrows();w++){
    int i = winPos[w];
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
    int i = winNeg[w];
    if (keepWin[i-1]==-1){
      probNegWin[i-1]=((1-2*proporNeg[w])/(1-proporNeg[w]))-errorRate;
      probPosWin[i-1]=errorRate;
    }
  }
  NumericVector probPosAl(posFragments.nrows()+1);  
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
  NumericVector probNegAl(negFragments.nrows()+1); 
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
  return List::create(
    _["Pos"] = listKeepPosReads,
    _["Neg"] = listKeepNegReads);
}
