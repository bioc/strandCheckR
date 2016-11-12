#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

//' @title  compute strand information of sliding window based on coverage
//'
//' @description non
//'
//' @param covPosLen the length the Rle object, which is the coverage comes from positive reads
//' @param covPosVal the value the Rle object, which is the coverage comes from positive reads
//' @param covNegLen the length the Rle object, which is the coverage comes from negative reads
//' @param covNegVal the value the Rle object, which is the coverage comes from negative reads
//' @param end the last base on the reference chromosome that the sliding window atteint
//' @param win the size of the sliding window
//' @param step the step of the sliding window
//' @param minCov if a window has the max coverage least than minCov, then it will be rejected
//' @param maxCov if a window has the max coverage greater than maxCov, then it will be kept
//' @param logitThreshold the logit of the threshold
//'
//' @return A list of two data frames Plus and Minus which respectively contains information of positive windows and negative windows: 'win' is the window number, and 'value' is the normalized estimated value to be tested
//' 
//' @examples
//' computeWinCov(c(10000,200,30,1,20,50),c(0,2,3,4,1,2),c(10020,300,20,1,10,15),c(0,1,3,5,1,3),10200,1000,100,2,50,0.8472979)
//' 
//' @export
// [[Rcpp::export]]

List computeWinCov(IntegerVector covPosLen,IntegerVector covPosVal,IntegerVector covNegLen,IntegerVector covNegVal,int end,int win,int step,int minCov,int maxCov,double logitThreshold){
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
    if (maxCovP>minCov || maxCovN>minCov){
      double estimate = (double)Plus/(Plus+Minus);
      double error = sqrt(1./(Plus+Minus)/estimate/(1-estimate));
      double lTestimate=log(estimate/(1-estimate));
      double value=(lTestimate - logitThreshold)/error;
      if (lTestimate<=0) value=-(lTestimate+logitThreshold)/error;
      if (Plus>Minus || (maxCovP>maxCov && maxCov>0)){
        if (Minus==0 || (maxCovP>maxCov && maxCov>0)) valueP.push_back(1e10);
        else valueP.push_back(value);
        windowP.push_back(c);
      }
      else if (Plus<=Minus || (maxCovN>maxCov & maxCov>0)){
        if (Plus==0 || (maxCovN>maxCov & maxCov>0)) valueM.push_back(1e10);
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



