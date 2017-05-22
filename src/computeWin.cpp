#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

//' @title  Compute strand information of sliding window 
//'
//' @description Compute the positive proportion and the value to be tested afterward to decide whether the window is kept or not (this value is calculated from the estimated proportion and error).
//' This method is called by the functions filterOne and filterMulti in the case that we don't need other plotting information from the windows.
//' 
//' @param covPosLen the run length of an Rle object which is the coverage comes from positive reads
//' @param covPosVal the run value of an Rle object which is the coverage comes from positive reads
//' @param covNegLen the run length of an Rle object which is the coverage comes from negative reads
//' @param covNegVal the run value of an Rle object which is the coverage comes from negative reads
//' @param end the last base on the reference chromosome that the sliding window atteints
//' @param readLength the average length of reads
//' @param win the size of the sliding window
//' @param step the step of the sliding window
//' @param minCov if a window has the max coverage smaller than minCov, then it will be rejected regardless its strand proportion.
//' @param maxCov if a window has the max coverage greater than maxCov, then it will be kept regardless its strand proportion. If maxCov=0 then it doesn't have any effect on selecting windows.
//'
//' @return A list of two data frames Plus and Minus which respectively contains information of positive windows and negative windows. 
//' Each data frame contains the information of window number, proportion of postive reads, and the value to be tested afterward to decide whether the window is kept or not (this value is calculated from the estimated proportion and error).
//' 
//' @seealso filterOne, filterMulti
//' 
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
//' maxCov <- 0
//' windows <- rnaCleanR::computeWin(runLength(covPos),runValue(covPos),runLength(covNeg),runValue(covNeg),readLength,len,win,step,minCov,maxCov)
//' 
//' @export
//' 
// [[Rcpp::export]]

List computeWin(IntegerVector covPosLen,IntegerVector covPosVal,IntegerVector covNegLen,IntegerVector covNegVal,double readLength,int end,int win,int step,int minCov,int maxCov){
  int start=0;
  int preP=0;
  int preM=0;
  int iP=0;
  int iM=0;
  std::vector<int> windowP;
  std::vector<int> windowM;
  std::vector<double> proporP;
  std::vector<double> proporM;
  std::vector<double> errorP;
  std::vector<double> errorM;
  int c=1;
  int maxCovP=0;
  int maxCovN=0;
  while (start<end){
    int Plus=increaseVal(covPosLen,covPosVal,iP,start,win,preP,maxCovP);//get the number of '+' in the current window
    int Minus=increaseVal(covNegLen,covNegVal,iM,start,win,preM,maxCovN);//get the number of '-' in the current window
    increase(covPosLen,iP,start,end,step,preP);//go to the next window of positive coverage
    increase(covNegLen,iM,start,end,step,preM);//go to the next window of negative coverage
    start+=step;
    if (Plus>0 || Minus>0){
      double propor = (double)Plus/(Plus+Minus);
      double max = maxCovP;
      if (propor<=0.5) max = maxCovN;
      if (max>minCov){
        double error = sqrt(readLength/(Plus+Minus)/propor/(1-propor));//The standard error of positive proportion upon logistic regression: 1/sqrt(n*p*(1-p)), here n = (Plus+Minus)/readLength
        if (Plus>Minus || (maxCovP>maxCov && maxCov>0)){
          windowP.push_back(c);
          proporP.push_back(propor);
          errorP.push_back(error);
        }
        if (Plus<=Minus || (maxCovN>maxCov & maxCov>0)){
          windowM.push_back(c);
          proporM.push_back(propor);
          errorM.push_back(error);
        }
      }
    }
    c++;
  }
  return List::create(
    _["Plus"] = DataFrame::create(_["win"]= windowP, _["propor"]= proporP, _["error"] = errorP),
      _["Minus"] = DataFrame::create(_["win"]= windowM, _["propor"]= proporM, _["error"] = errorM)
  );
}



