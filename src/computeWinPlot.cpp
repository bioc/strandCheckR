#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

//' @title  Compute strand information of sliding window (plot version)
//'
//' @description Compute the positive proportion, sum of reads, max coverage and the group of each window. Windows are grouped based on their maximum coverage: by default definition, groups spead from 1 to 8, which correspond to the max coverage respectively in the range "0-10","10-20","20-50","50-100","100-200","200-500","500-1000",">1000"
//' This method is used in the getPlot function when we only need the information to plot, and do not need to filter the reads afterward.
//' 
//' @param covPosLen the run length of an Rle object which is the coverage comes from positive reads
//' @param covPosVal the run value of an Rle object which is the coverage comes from positive reads
//' @param covNegLen the run length of an Rle object which is the coverage comes from negative reads
//' @param covNegVal the run value of an Rle object which is the coverage comes from negative reads
//' @param end the last base on the reference chromosome that the sliding window atteint
//' @param readLength the average length of reads
//' @param win the size of the sliding window
//' @param step the step of the sliding window
//' @param minCov if a window has the max coverage least than minCov, then it will not be counted
//'
//' @return A list of two data frames Plus and Minus which respectively contains information of positive windows and negative windows: 'win' is the window number, and 'value' is the normalized estimated value to be tested
//' Each data frame contains contain the information of proportion of postive reads, the max coverage and the group of max coverage
//' @examples
//' bamfilein <- system.file("data","s1.chr1.bam",package = "rnaCleanR")
//' alignment <- GenomicAlignments::readGAlignments(bamfilein) 
//' alignmentInChr1 <- alignment[seqnames(alignment)=="1"] 
//' covPos <- alignmentInChr1[strand(alignment)=="+"] %>% GenomicAlignments::coverage() 
//' covNeg <- alignmentInChr1[strand(alignment)=="-"] %>% GenomicAlignments::coverage() 
//' len <- length(covChr)
//' readLength <- 100
//' win <- 1000
//' step <- 100
//' minCov <- 0
//' windows <- computeWinPlot(runLength(covPos),runValue(covPos),runLength(covNeg),runValue(covNeg),readLength,len,win,step,minCov)
//' 
//' 
//' @export
//' 
// [[Rcpp::export]]

List computeWinPlot(IntegerVector covPosLen,IntegerVector covPosVal,IntegerVector covNegLen,IntegerVector covNegVal,double readLength,int end,int win,int step,int minCov){
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
        window.push_back(c);
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
  return DataFrame::create(_["win"]=window,_["propor"]=pro, _["sum"]= sum, _ ["max"]= maxC, _ ["group"] = group);
}



