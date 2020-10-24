
#include <Rcpp.h>
#include <ANN/ANN.h>     // ANN library header
#include "WANN.h"

using namespace Rcpp;

// [[Rcpp::export]]
List nn2_cpp(NumericMatrix data, NumericMatrix query, const int k, const double eps = 0.0, const double radius = NA_REAL) {

  WANN tree = WANN(data);
  // NB this method will do a regular search when radius=0 or NA
  return tree.query_FR(query, k, radius, eps);
}
