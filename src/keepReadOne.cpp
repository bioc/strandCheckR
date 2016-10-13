#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List keepReadOne(IntegerVector startPosAl,IntegerVector endPosAl,IntegerVector groupNamePos,IntegerVector startNegAl,IntegerVector endNegAl,IntegerVector groupNameNeg,IntegerVector keepWinPos,IntegerVector keepWinNeg,int end,int win,int step,double limit){
  int minP = keepWinPos[0];
  int maxP = keepWinPos[keepWinPos.size()-1];
  bool* keepPos = new bool[maxP+1];
  for (int i=0;i<=maxP;i++) keepPos[i]=false;
  for (int i=0;i<keepWinPos.size();i++){
    keepPos[keepWinPos[i]]=true;
  }
  int minM = keepWinNeg[0];
  int maxM = keepWinNeg[keepWinNeg.size()-1];
  bool* keepNeg = new bool[maxM+1];
  for (int i=0;i<=maxM;i++) keepNeg[i]=false;
  for (int i=0;i<keepWinNeg.size();i++){
    keepNeg[keepWinNeg[i]]=true;
  }
  std::vector<int> keepReadPos;
  std::vector<int> keepReadNeg;
  int countP=0;
  int countM=0;
  for (int a=0;a<startPosAl.size();a++){
    int s = startPosAl[a];
    int e = endPosAl[a];
    int lim = (e-s+1)*limit;
    s = s + lim;
    e = e - lim;
    int wS=0;
    if (s>win){
      wS=ceil((s-win-1)/(double)step);//first windows that contain fragment i
    }
    int wE=maxP;
    if (e<(1+step*(maxP))){
      wE=floor(e/(double)step);//last window that contain fragment i
    }
    bool bl=false;
    if (wE>=minP && wS<=maxP){
      int i=wS;
      if (i<minP) i=minP;
      int j=wE;
      if (j>maxP) j=maxP;
      bl = keepPos[i];
      i++;
      while (!bl && i<=j){
        bl = keepPos[i];
        i++;
      }
      if (bl){ keepReadPos.push_back(groupNamePos[a]);}
    }
  }
  for (int a=0;a<startNegAl.size();a++){
    int s = startNegAl[a];
    int e = endNegAl[a];
    int limit = (e-s+1)*25/100;
    s = s + limit;
    e = e - limit;
    int wS=0;
    if (s>win){
      wS=ceil((s-win-1)/(double)step);//first windows that contain fragment i
    }
    int wE=maxM;
    if (e<(1+step*(maxM))){
      wE=floor(e/(double)step);//last window that contain fragment i
    }
    bool bl=false;
    if (wE>=minM && wS<=maxM){
      int i=wS;
      if (i<minM) i=minM;
      int j=wE;
      if (j>maxM) j=maxM;
      bl = keepNeg[i];
      i++;
      while (!bl && i<=j){
        bl = keepNeg[i];
        i++;
      }
      if (bl){ keepReadNeg.push_back(groupNameNeg[a]);}
    }
  }
  delete[] keepPos;
  delete[] keepNeg; 
  return List::create(
      _["Pos"] = keepReadPos,
     _["Neg"] = keepReadNeg
  );
}
