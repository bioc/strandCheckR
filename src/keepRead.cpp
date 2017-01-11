#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List keepRead(IntegerVector startPosAl,IntegerVector endPosAl,IntegerVector groupPos,IntegerVector startNegAl,IntegerVector endNegAl,IntegerVector groupNeg,NumericVector proporWinPos,NumericVector proporWinNeg,IntegerVector keepWinPos,IntegerVector keepWinNeg,int maxWin,int win,int step,double errorRate){
  int minP = startPosAl[0];
  int minN = startNegAl[0];
  int* keepWin = new int[maxWin];
  for (int i=0;i<maxWin;i++) keepWin[i]=0;
  for (int i=0;i<keepWinPos.size();i++){
    keepWin[keepWinPos[i]-1]=1;
  }
  for (int i=0;i<keepWinNeg.size();i++){
    if (keepWin[keepWinNeg[i]-1]==1){
      keepWin[keepWinNeg[i]-1]=2;
    }
    else keepWin[keepWinNeg[i]-1]=-1;
  }
  NumericVector probPosWin(maxWin);
  NumericVector probNegWin(maxWin);
  for (int i=0;i<maxWin;i++){
    probPosWin[i] = errorRate;
    probNegWin[i] = errorRate;
  }
  for (int w=0;w<keepWinPos.length();w++){
    int i = keepWinPos[w];
    if (keepWin[i-1]==1){
      probPosWin[i-1]=(2-1/proporWinPos[w]);
    }
    else if (keepWin[i-1]==2){
      probPosWin[i-1]=(1-errorRate);
      probNegWin[i-1]=(1-errorRate);
    }
  }
  for (int w=0;w<keepWinNeg.length();w++){
    int i = keepWinNeg[w];
    if (keepWin[i-1]==-1){
      probNegWin[i-1]=((1-2*proporWinNeg[w])/(1-proporWinNeg[w]));
    }
  }
  NumericVector probPosAl(startPosAl.size()+1);  
  for (int i=0;i<startPosAl.size();i++) probPosAl[i]=0;
  for (int a=0;a<startPosAl.size();a++){
    int s = startPosAl[a];
    int e = endPosAl[a];
    int x1 = 1+ceil((double)(s-win)/step);//the first windows that the start position of the read is in
    int x2 = 1+floor((double)(e-1)/step);//the last windows that the end position of the read is in
    if (x1<0) x1=0;
    double prob = 0;
    for (int k=x1;k<=x2;k++){ 
      prob=std::max(prob,probPosWin[k]);
    }
    probPosAl[groupPos[a]]=std::max(probPosAl[groupPos[a]],prob);
  }
  NumericVector probNegAl(startNegAl.size()+1); 
  for (int i=0;i<startNegAl.size();i++) probNegAl[i]=0;
  for (int a=0;a<startNegAl.size();a++){
    int s = startNegAl[a];
    int e = endNegAl[a];
    int x1 = 1+ceil((double)(s-win)/step);//the first windows that the start position of the read is in
    int x2 = 1+floor((double)(e-1)/step);//the last windows that the end position of the read is in
    if (x1<0) x1=0;
    double prob = 0;
    for (int k=x1;k<=x2;k++){ 
      prob=std::max(prob,probNegWin[k]);
    }
    probNegAl[groupNeg[a]]=std::max(probNegAl[groupNeg[a]],prob);
  }
  
  std::default_random_engine generator;
  std::vector<int> listKeepPosReads;
  std::vector<int> listKeepNegReads;
  for (int i=0;i<startPosAl.size();i++){
    std::binomial_distribution<int> distribution(1,probPosAl[groupPos[i]]);
    if (distribution(generator)){
      listKeepPosReads.push_back(groupPos[i]);
    }
  }
  for (int i=0;i<startNegAl.size();i++){
    std::binomial_distribution<int> distribution(1,probNegAl[groupNeg[i]]);
    if (distribution(generator)){
      listKeepNegReads.push_back(groupNeg[i]);
    }
  }
  return List::create(
    _["Pos"] = listKeepPosReads,
    _["Neg"] = listKeepNegReads);
}
