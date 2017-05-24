#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

//' @title  Compute strand information of sliding window (plot version)
//'
//' @description Compute the positive proportion, sum of reads, max coverage and the group of each window. Windows are grouped based on their maximum coverage: by default definition, groups spead from 1 to 4, which correspond to the max coverage respectively in the range "0-10","10-100","100-1000",">1000"
//' This method is used in the getPlot function when we only need the information to plot, and do not need to filter the reads afterward.
//' 
//' @param covPosLen the run length of an Rle object which is the coverage comes from positive reads
//' @param covPosVal the run value of an Rle object which is the coverage comes from positive reads
//' @param covNegLen the run length of an Rle object which is the coverage comes from negative reads
//' @param covNegVal the run value of an Rle object which is the coverage comes from negative reads
//' @param end the last base on the reference chromosome that the sliding window atteints
//' @param readLength the average length of reads
//' @param win the size of the sliding window
//' @param step the step of the sliding window
//' @param minCov if a window has the max coverage least than minCov, then it will not be counted
//'
//' @return A data frame which contains the information of all windows: Starting positive, Number of Positive/Negative Reads, the Maximum Coverage.
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
//' windows <- rnaCleanR::computeWinPlot(runLength(covPos),runValue(covPos),runLength(covNeg),runValue(covNeg),readLength,len,win,step,minCov)
//' 
//' 
//' @export
//' 
// [[Rcpp::export]]

List computeWinInfo(IntegerVector covPosLen,IntegerVector covPosVal,IntegerVector covNegLen,IntegerVector covNegVal,double readLength,int end,int win,int step,int minCov){
  int start=0;
  int preP=0;
  int preM=0;
  int iP=0;
  int iM=0;
  
  std::vector<int> Start;
  std::vector<double> PositiveReads;
  std::vector<double> NegativeReads;
  std::vector<int> MaxC;
  int maxCovP=0;
  int maxCovN=0;
  while (start<end){
    int Plus=increaseVal(covPosLen,covPosVal,iP,start,win,preP,maxCovP);
    Plus = round(Plus/readLength);
    int Minus=increaseVal(covNegLen,covNegVal,iM,start,win,preM,maxCovN);
    Minus = round(Minus/readLength);
    increase(covPosLen,iP,start,end,step,preP);
    increase(covNegLen,iM,start,end,step,preM);
    if (Plus>0 || Minus>0){
      double max = maxCovP+maxCovN;
      if (max>minCov){
        Start.push_back(start+1);
        PositiveReads.push_back(Plus);
        NegativeReads.push_back(Minus);
        MaxC.push_back(max);
      }
    }
    start+=step;
  }
  return DataFrame::create(_["Start"]=Start,_["NbPositiveReads"]=PositiveReads, _["NbNegativeReads"]= NegativeReads, _ ["MaxCoverage"]= MaxC);
}



