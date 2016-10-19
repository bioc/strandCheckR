#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
List computeWinCount0(IntegerVector startPos,IntegerVector endPos,IntegerVector startNeg,IntegerVector endNeg,int end,int win,int step,double limit,double logitThreshold){
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
  for (int i=0;i<nbWin;i++){
    windowP[i]=0;
    windowM[i]=0;
  }
  int i=0;
  while (i<startPos.size() && startPos[i]<end){
    int s = startPos[i];
    int e = endPos[i];
    int lim = (int)(e-s+1)*limit;
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
      }
    }
    i++;
  }
  i=0;
  while (i<startNeg.size() && startNeg[i]<end){
    int s = startNeg[i];
    int e = endNeg[i];
    int lim = (int)(e-s+1)*limit;
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
      }   
    }
    i++;
  }
  std::vector<double> valueP;
  std::vector<int> winP;
  std::vector<double> valueM;
  std::vector<int> winM;
    for (int i=0;i<nbWin;i++){
      int Plus=windowP[i];
      int Minus=windowM[i];
      if ((Plus>0 || Minus>0)){
        double estimate = (double)Plus/(Plus+Minus);
        double error = sqrt(1./(Plus+Minus)/estimate/(1-estimate));
        double lTestimate=log(estimate/(1-estimate));
        double value=(lTestimate - logitThreshold)/error;
        if (lTestimate<=0) value=-(lTestimate+logitThreshold)/error;
        if (Plus>Minus){
          if (Minus==0) valueP.push_back(1e10);
          else valueP.push_back(value);
          winP.push_back(i+1);
        }
        else{
          if (Plus==0) valueM.push_back(1e10);
          else valueM.push_back(value);
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


