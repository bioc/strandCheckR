#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
List computeWinCov(IntegerVector covPosLen,IntegerVector covPosVal,IntegerVector covNegLen,IntegerVector covNegVal,int end,int win,int step,double logitThreshold,int minR,int maxR){
  int start=0;
  int preP=0;
  int preM=0;
  int iP=0;
  int iM=0;
  std::vector<int> windowP;
  std::vector<int> windowM;
  std::vector<double> valueP;
  std::vector<double> valueM;
  int c=1;
  int maxCovP=0;
  int maxCovN=0;
  while (start<end){
    int Plus=increaseVal(covPosLen,covPosVal,iP,start,win,preP,maxCovP);
    int Minus=increaseVal(covNegLen,covNegVal,iM,start,win,preM,maxCovN);
    increase(covPosLen,iP,start,end,step,preP);
    increase(covNegLen,iM,start,end,step,preM);
    start+=step;
    if (maxCovP>minR || maxCovN>minR){
      double estimate = (double)Plus/(Plus+Minus);
      double error = sqrt(1./(Plus+Minus)/estimate/(1-estimate));
      double lTestimate=log(estimate/(1-estimate));
      double value=(lTestimate - logitThreshold)/error;
      if (lTestimate<=0) value=-(lTestimate+logitThreshold)/error;
      if (Plus>Minus || maxCovP>maxR){
        if (Minus==0 || maxCovP>maxR) valueP.push_back(1e10);
        else valueP.push_back(value);
        windowP.push_back(c);
      }
      else if (Plus<=Minus || maxCovN<maxR){
        if (Plus==0 || maxCovN<maxR) valueM.push_back(1e10);
        else valueM.push_back(value);
        windowM.push_back(c);
      }
    }
    c++;
  }
  return List::create(
    _["Plus"] = DataFrame::create(_["win"]= windowP, _["value"]= valueP),
    _["Minus"] = DataFrame::create(_["win"]= windowM, _["value"]= valueM)
  );
}



