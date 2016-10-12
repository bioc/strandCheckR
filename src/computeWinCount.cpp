#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
List computeWinCount(IntegerVector startPos,IntegerVector endPos,IntegerVector startNeg,IntegerVector endNeg,int end,int win,int step,double logitThreshold,int minR){
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
    int limit = (e-s+1)*25/100;
    s = s + limit;
    e = e - limit;
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
    while (currentWin<wS|| (i==startPos.size()-1 && currentWin<=wE)){//current window does not contain fragment i
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
    int limit = (e-s+1)*25/100;
    s = s + limit;
    e = e - limit;
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
    while (currentWin<wS || (i==startNeg.size()-1 && currentWin<=wE)){//current window does not contain fragment i
      if (windowM[currentWin]>0){
        maxM[currentWin]=getMaxCov(startNeg,endNeg,fragsInWin[currentWin],currentWin,win,step);
        fragsInWin[currentWin].clear();
      }
      currentWin++;
    }
    i++;
}
  std::vector<double> valueP;
  std::vector<int> winP;
  std::vector<int> maxCovP;
  std::vector<int> countP;
  std::vector<double> perP;
  std::vector<double> valueM;
  std::vector<int> winM;
  std::vector<int> maxCovM;
  std::vector<int> countM;
  std::vector<double> perM;
  for (int i=0;i<nbWin;i++){
    int Plus=windowP[i];
    int Minus=windowM[i];
    if ((maxP[i]>minR || maxM[i]>minR)){
      double estimate = (double)Plus/(Plus+Minus);
      double error = sqrt(1./(Plus+Minus)/estimate/(1-estimate));
      double lTestimate=log(estimate/(1-estimate));
      double value=(lTestimate - logitThreshold)/error;
      if (lTestimate<=0) value=-(lTestimate+logitThreshold)/error;
      if (Plus>Minus){
        if (Minus==0) valueP.push_back(1e10);
        else valueP.push_back(value);
        winP.push_back(i+1);
        countP.push_back(Plus);
        perP.push_back((double)Plus/(Plus+Minus));
        maxCovP.push_back(maxP[i]);
      }
      else{
        if (Plus==0) valueM.push_back(1e10);
        else valueM.push_back(value);
        winM.push_back(i+1);
        countM.push_back(Minus);
        perM.push_back((double)Minus/(Plus+Minus));
        maxCovM.push_back(maxM[i]);
      }
    }
  }
  return List::create(
    _["Plus"] = DataFrame::create(_["win"]= winP,//_["maxCov"]=maxCovP, _["count"]=countP, _["percent"]=perP, 
                              _["value"]= valueP),
      _["Minus"] =DataFrame::create(_["win"]= winM,//_["maxCov"]=maxCovM, _["count"]=countM, _["percent"]=perM, 
                               _["value"]= valueM)
  );
}


