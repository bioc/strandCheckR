#include <Rcpp.h>
using namespace Rcpp;

//' @title  Get the indices of positive and negative reads in a given alignment
//'
//' @description non
//'
//' @param strand a string vector contain the "+" or "_"
//'
//' @return A list of two vectors which contains the order of each positive, negative strands within the input strand vector
//' 
//' @examples
//' strand <- rep("+",1000)
//' negInd <- sample(1:1000,400,replace = FALSE)
//' strand[negInd] <- "_"
//' index <- getIndex(strand)
//' 
//' 
//' @export
//' 
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
