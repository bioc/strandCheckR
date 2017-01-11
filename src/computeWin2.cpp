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
//' @param readLength the average length of reads
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
//' computeWin2(c(10000,200,30,1,20,50),c(0,2,3,4,1,2),c(10020,300,20,1,10,15),c(0,1,3,5,1,3),100,10200,1000,100,0)
//' 
//' @export
// [[Rcpp::export]]

List computeWin2(IntegerVector covPosLen,IntegerVector covPosVal,IntegerVector covNegLen,IntegerVector covNegVal,double readLength,int end,int win,int step,int minCov){
  int start=0;
  int preP=0;
  int preM=0;
  int iP=0;
  int iM=0;
  
  std::vector<int> window;
  std::vector<double> pro;
  std::vector<double> sum;
  std::vector<int> maxC;
  std::vector<int> group;
  int c=1;
  int maxCovP=0;
  int maxCovN=0;
  while (start<end){
    int Plus=increaseVal(covPosLen,covPosVal,iP,start,win,preP,maxCovP);
    int Minus=increaseVal(covNegLen,covNegVal,iM,start,win,preM,maxCovN);
    increase(covPosLen,iP,start,end,step,preP);
    increase(covNegLen,iM,start,end,step,preM);
    start+=step;
    if (Plus>0 || Minus>0){
      double estimate = (double)Plus/(Plus+Minus);
      double max = maxCovP;
      if (estimate<=0.5) max = maxCovN;
      if (max>minCov){
        pro.push_back(estimate);
        sum.push_back((Plus+Minus)/(double)readLength);
        maxC.push_back(max);
        if (max>1000) group.push_back(8);
        else if (max>500) group.push_back(7);
        else if (max>200) group.push_back(6);
        else if (max>100) group.push_back(5);
        else if (max>50) group.push_back(4);
        else if (max>20) group.push_back(3);
        else if (max>10) group.push_back(2);
        else if (max>0) group.push_back(1);
      }
      
    }
    c++;
  }
  return DataFrame::create(_["propor"]=pro, _["sum"]= sum, _ ["max"]= maxC, _ ["group"] = group);
}



