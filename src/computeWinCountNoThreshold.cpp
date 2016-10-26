#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
List computeWinCountNoThreshold(IntegerVector startPos,IntegerVector endPos,IntegerVector startNeg,IntegerVector endNeg,int end,int win,int step,int minR,int maxR,double limit){
  int start=0;
  int preP=0;
  int preM=0;
  int iP=0;
  int iM=0;
  int nbWin = 0;
  if (end<=win) nbWin=1;
  else nbWin=1+ceil((end-win-1)/(double)step);
  IntegerVector windowP(nbWin);
  IntegerVector windowM(nbWin);
  IntegerVector maxP(nbWin);
  IntegerVector maxM(nbWin);
  for (int i=0;i<nbWin;i++){
    windowP[i]=0;
    windowM[i]=0;
    maxP[i]=0;
    maxM[i]=0;
  }
  std::vector<int>* fragsInWin = new std::vector<int>[nbWin];
  int currentWin=0;
  int i=0;
  while (i<startPos.size() && startPos[i]<end){
    int s = startPos[i];
    int e = endPos[i];
    int lim = (e-s+1)*limit;
    s = s + lim;
    e = e - lim;
    int wS=0;
    if (s>win){
      wS=ceil((s-win-1)/(double)step);//first windows that contain fragment i
    }
    int wE=nbWin-1;
    if (e<(1+step*(nbWin-1))){
      wE=floor(e/(double)step);//last window that contain fragment i
    }
    if (wS<=wE) {
      for (int j=wS;j<=wE;j++){ 
        windowP[j]++;
        fragsInWin[j].push_back(i);
      }
    }
    ///Compute max coverage
    while (currentWin<wS|| (i==startPos.size()-1 && currentWin<nbWin)){//current window does not contain fragment i
      if (windowP[currentWin]>0){
        maxP[currentWin]=getMaxCov(startPos,endPos,fragsInWin[currentWin],currentWin,win,step);
        fragsInWin[currentWin].clear();
      }
      currentWin++;
    }
    i++;
  }
  currentWin=0;
  i=0;
  while (i<startNeg.size() && startNeg[i]<end){
    int s = startNeg[i];
    int e = endNeg[i];
    int lim = (e-s+1)*limit;
    s = s + lim;
    e = e - lim;
    int wS=0;
    if (s>win){
      wS=ceil((s-win-1)/(double)step);
    }
    int wE=nbWin-1;
    if (e<(1+step*(nbWin-1))){
      wE=floor(e/(double)step);
    }
    if (wS<=wE){
      for (int j=wS;j<=wE;j++){
        windowM[j]++;
        fragsInWin[j].push_back(i);
      }   
    }
    //Compute max coverage
    while (currentWin<wS || (i==startNeg.size()-1 && currentWin<nbWin)){//current window does not contain fragment i
      if (windowM[currentWin]>0){
        maxM[currentWin]=getMaxCov(startNeg,endNeg,fragsInWin[currentWin],currentWin,win,step);
        fragsInWin[currentWin].clear();
      }
      currentWin++;
    }
    i++;
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
        cP++;
      }
      else if (maxP[i]>50) cnP++;
    }
    else{
      if (maxM[i]>=5 && maxM[i]<=50){
        sM+=(double)Minus/(Plus+Minus);
        cM++;
      }
      else if (maxM[i]>50) cnM++;
    }
  }
  double thresholdP=(sP/cP);
  if (cP/(double)(cP+cnP)>=0.9) std::cout<<"Suggested threshold for positive window is "<<thresholdP<<std::endl;
  else {
    thresholdP=0.7;
    std::cout<<"Can not estimate threshold for positive window, take the default threshold "<<thresholdP<<std::endl;
  }
  double thresholdM=(sM/cM);
  if (cM/(double)(cM+cnM)>=0.9) std::cout<<"Suggested threshold for negative window is  "<<thresholdM<<std::endl;
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
    if ((maxP[i]>minR || maxM[i]>minR)){
      double estimate = (double)Plus/(Plus+Minus);
      double error = sqrt(1./(Plus+Minus)/estimate/(1-estimate));
      double lTestimate=log(estimate/(1-estimate));
      if (Plus>Minus || (maxP[i]>maxR && maxR>0)){
        if (Minus==0 || (maxP[i]>maxR && maxR>0)) valueP.push_back(1e10);//always keep these windows
        else {
          double value=(lTestimate - logitThresholdP)/error;
          if (lTestimate<=0) value=-(lTestimate+logitThresholdP)/error;
          valueP.push_back(value);
        }
        winP.push_back(i+1);
      }
      else if (Plus<=Minus || (maxM[i]>maxR && maxR>0)){
        if (Plus==0 || (maxM[i]>maxR && maxR>0)) valueM.push_back(1e10);//always keep these windows
        else{
          double value=(lTestimate - logitThresholdM)/error;
          if (lTestimate<=0) value=-(lTestimate+logitThresholdM)/error;
          valueM.push_back(value);
        } 
        winM.push_back(i+1);
      }
    }
  }
  return List::create(
    _["Plus"] = DataFrame::create(_["win"]= winP,
                              _["value"]= valueP),
    _["Minus"] =DataFrame::create(_["win"]= winM,
                              _["value"]= valueM)
  );
}


