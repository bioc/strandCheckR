#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List getIndex(CharacterVector strand){
  std::vector<int> indexPos;
  std::vector<int> indexNeg;
  for (int i=0;i<strand.size();i++){
      if (strand[i]=="+"){
        indexPos.push_back(i+1);
      } 
      else{
        indexNeg.push_back(i+1);
      }
  }
  return List::create(
    _["Pos"] = indexPos,
    _["Neg"] = indexNeg
  );
}
