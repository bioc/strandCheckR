#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

void increase(IntegerVector len,int& ind, int start, int end, int inc, int& pre){
  int next = start;
  while (ind<len.size() && next < start + inc){
    int s = next;
    next = s+len[ind]-pre;
    if (next>start+inc){
      pre+=(start+inc)-s;
      next=start+inc;
    }
    else{
      pre=0;
      ind++;
    }
  }
  start=next;
}

int increaseVal(IntegerVector len,IntegerVector val,int ind, int start, int inc, int pre,int& maxCov){
  int next = start;
  int sum = 0;
  maxCov = 0;
  while (ind<len.size() && next < start + inc){
    int s = next;
    next = s+len[ind]-pre;
    int v = val[ind];
    if (v>maxCov) maxCov=v;
    if (next>start+inc){
      pre+=(start+inc)-s;
      next=start+inc;
    }
    else{
      pre=0;
      ind++;
    }
    sum+=v*(next-s);
  }
  return sum;
}


int getMaxCov(IntegerVector start,IntegerVector end,std::vector<int> fragsInWin,int indiceW,int win,int step){
  std::vector<int>::iterator iter = fragsInWin.begin();
  int dem = indiceW*step+1;
  NumericVector cov(win);
  for (int i =0;i<win;i++) cov(i)=0;
  while (iter!=fragsInWin.end()){
    for (int i=start(*iter)-dem;i<=end(*iter)-dem;i++){
      if (i>=0 && i<win) {
        cov(i)++;
      }
    }
    iter++;
  }
  int max=0;
  for (int i=0;i<win;i++) if (cov(i)>max) max=cov(i);
  return max;
}
