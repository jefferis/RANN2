
#include <Rcpp.h>
#include "ANN.h"     // ANN library header

using namespace Rcpp;

// [[Rcpp::export]]
List nn2(NumericMatrix data, NumericMatrix query, int k) {

	List z = List::create( data, query ) ;
	return z ;
}
