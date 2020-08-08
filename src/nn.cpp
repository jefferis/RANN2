
#include <Rcpp.h>
#include <ANN/ANN.h>     // ANN library header
#include "WANN.h"

using namespace Rcpp;

// [[Rcpp::export]]
List nn2_cpp(NumericMatrix data, NumericMatrix query, const int k, const double eps = 0.0, const double radius = NA_REAL) {

  WANN tree = WANN(data);

  if(R_IsNA(radius)) {
    return tree.query(query, k, eps);
  } else {
    return tree.query_FR(query, k, radius, eps);
  }
}
