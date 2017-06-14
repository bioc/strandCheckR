#ifndef PKG_UTILS_H
#define PKG_UTILS_H

using namespace Rcpp;
void increase(IntegerVector len,int& ind, int start, int end, int inc, int& pre);
int increaseVal(IntegerVector len,IntegerVector val,int ind, int start, int inc, int pre,int& maxCov);
int getMaxCov(IntegerVector start,IntegerVector end,std::vector<int> fragsInWin,int indiceW,int win,int step);
int increaseValR(IntegerVector len,DoubleVector val,int ind,int start, int end, int pre);
void increaseR(IntegerVector len,int& ind,int start, int end, int& pre);
#endif
