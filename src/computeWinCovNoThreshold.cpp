#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
List computeWinCovNoThreshold(IntegerVector covPosLen,IntegerVector covPosVal,IntegerVector covNegLen,IntegerVector covNegVal,int end,int win,int step,int minR,int maxR){
  int start=0;
  int preP=0;
  int preM=0;
  int iP=0;
  int iM=0;
  int nbWin = 0;
  if (end<=win) nbWin=1;
  else nbWin=1+ceil((end-1)/(double)step);
  IntegerVector windowP(nbWin);
  IntegerVector windowM(nbWin);
  IntegerVector maxP(nbWin);
  IntegerVector maxM(nbWin);
  int maxCovP=0;
  int maxCovN=0;
  int c=0;
  while (start<end){
    int Plus=increaseVal(covPosLen,covPosVal,iP,start,win,preP,maxCovP);
    int Minus=increaseVal(covNegLen,covNegVal,iM,start,win,preM,maxCovN);
    increase(covPosLen,iP,start,end,step,preP);
    increase(covNegLen,iM,start,end,step,preM);
    start+=step;
    windowP[c] = Plus;
    windowM[c] = Minus;
    maxP[c] = maxCovP;
    maxM[c] = maxCovN;
    c++;
  }
  double sP=0;
  int cP=0;
  int cnP=0;
  double sM=0;
  int cM=0;
  int cnM=0;
  for (int i=0;i<nbWin;i++){
    int Plus=windowP[i];
    int Minus=windowM[i];
    if (Plus>Minus){
      if (maxP[i]>=5 && maxP[i]<=50){
        sP+=(double)Plus/(Plus+Minus);
        cP+=1;
      }
      else if (maxP[i]>50) cnP++;
    }
    else{
      if (maxM[i]>=5 && maxM[i]<=50){
        sM+=(double)Minus/(Plus+Minus);
        cM+=1;
      }
      else if (maxM[i]>50) cnM++;
    }
  }
  double thresholdP=(sP/cP);
  if (cP/(double)(cP+cnP)>=0.95) std::cout<<"Suggested threshold for positive window is "<<thresholdP<<std::endl;
  else {
    thresholdP=0.7;
    std::cout<<"Can not estimate threshold for positive window, take the default threshold "<<thresholdP<<std::endl;
  }
  double thresholdM=(sM/cM);
  if (cM/(double)(cM+cnM)>=0.95) std::cout<<"Suggested threshold for negative window is "<<thresholdM<<std::endl;
  else{
    thresholdM=0.7;
    std::cout<<"Can not estimate threshold for negative window, take the default threshold "<<thresholdP<<std::endl;
  }
  double logitThresholdP=log(thresholdP/(1-thresholdP));
  double logitThresholdM=log(thresholdM/(1-thresholdM));
  std::vector<double> valueP;
  std::vector<int> winP;
  std::vector<double> valueM;
  std::vector<int> winM;
  for (int i=0;i<nbWin;i++){
    int Plus=windowP[i];
    int Minus=windowM[i];
    maxCovP=maxP[i];
    maxCovN=maxM[i];
    if (maxCovP>minR || maxCovN>minR){
      double estimate = (double)Plus/(Plus+Minus);
      double error = sqrt(1./(Plus+Minus)/estimate/(1-estimate));
      double lTestimate=log(estimate/(1-estimate));
      if (Plus>Minus || maxCovP>maxR){
        double value=(lTestimate - logitThresholdP)/error;
        if (lTestimate<=0) value=-(lTestimate+logitThresholdP)/error;
        if (Minus==0 || (maxCovP>maxR && maxR>0)) valueP.push_back(1e10);
        else valueP.push_back(value);
        winP.push_back(i+1);
      }
      else if (Plus<=Minus || maxCovN>maxR){
        double value=(lTestimate - logitThresholdM)/error;
        if (lTestimate<=0) value=-(lTestimate+logitThresholdM)/error;
        if (Plus==0 || (maxCovN>maxR && maxR>0)) valueM.push_back(1e10);
        else valueM.push_back(value);
        winM.push_back(i+1);
      }
    }
  }
  return List::create(
    _["Plus"] = DataFrame::create(_["win"]= winP, _["value"]= valueP),
    _["Minus"] = DataFrame::create(_["win"]= winM, _["value"]= valueM)
  );
}



