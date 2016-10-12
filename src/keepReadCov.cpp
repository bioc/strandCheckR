#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector keepReadCov(IntegerVector startAl,IntegerVector endAl,CharacterVector strand,IntegerVector keepWinPos,IntegerVector keepWinNeg,int end,int win,int step){
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
  std::vector<int> keepReads;
  for (int a=0;a<startAl.size();a++){
    int s = startAl[a];
    int e = endAl[a];
    int x1 = 1+ceil((double)(s-win)/step);//the first windows that the start position of the read is in
    int x2 = 1+floor((double)(s-1)/step);//the last windows that the start position of the read is in
    int y1 = 1+ceil((double)(e-win)/step);//the first windows that the end position of the read is in
    int y2 = 1+floor((double)(e-1)/step);//the last windows that the end position of the read is in
    bool bl=false;
    if (strand[a]=="+"){
      if (x2>=minP && x1<=maxP){
        int i=x1;
        if (i<minP) i=minP;
        int j=x2;
        if (j>maxP) j=maxP;
        bl = keepPos[i];
        while (!bl && i<=j){
          i++;
          bl = keepPos[i];
        }
        if (bl){ keepReads.push_back(a+1);}
      }
      if (!bl){
        if (y1<=maxP && y2>=minP){
          int i=y1;
          if (i<minP) i=minP;
          int j=y2;
          if (j>maxP) j=maxP;
          bool bl = keepPos[i];
          while (!bl && i<=j){
            i++;
            bl = keepPos[i];
          }
          if (bl){ keepReads.push_back(a+1);}
        }
      }
    }
    else{
      if (x2>=minM && x1<=maxM){
        int i=x1;
        if (i<minM) i=minM;
        int j=x2;
        if (j>maxM) j=maxM;
        bl = keepNeg[i];
        while (!bl && i<=j){
          i++;
          bl = keepNeg[i];
        }
        if (bl){ keepReads.push_back(a+1);}
      }
      if (!bl){
        if (y1<=maxM && y2>=minM){
          int i=y1;
          if (i<minM) i=minM;
          int j=y2;
          if (j>maxM) j=maxM;
          bool bl = keepNeg[i];
          while (!bl && i<=j){
            i++;
            bl = keepNeg[i];
          }
          if (bl){ keepReads.push_back(a+1);}
        }
      }
    }
  }
  delete[] keepPos;
  delete[] keepNeg;
  return wrap(keepReads);
}
