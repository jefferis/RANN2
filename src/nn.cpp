
#include <Rcpp.h>
#include <ANN/ANN.h>     // ANN library header
#include "WANN.h"

using namespace Rcpp;

// [[Rcpp::export]]
List nn2_cpp(NumericMatrix data, NumericMatrix query, const int k, const double eps = 0.0, const double radius = NA_REAL) {

  WANN tree = WANN(data);

  if(!R_IsNA(radius) && radius < 0.0){
    // this will return a ragged list of vectors allowing searches with highly variable numbers of neighbors
    return tree.query_FR_ragged(query, k, abs(radius), eps);
  }else{
    // NB this method will do a regular search when radius=0 or NA
    return tree.query_FR(query, k, radius, eps);
  }
}
