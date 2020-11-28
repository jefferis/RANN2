#ifndef WANN_H
#define WANN_H

#include <Rcpp.h>
#include <ANN/ANN.h>     // ANN library header
using namespace Rcpp;

class WANN {
private:
  ANNpointArray data_pts;
  ANNkd_tree  *tree;
  int d;
  int n;

public:
  WANN(NumericMatrix data, bool buildtree=true);

  ~WANN() {
    annDeallocPts(data_pts);
    delete_tree();
  }

  void build_tree() {
    if(tree==0) {
      tree = new ANNkd_tree( data_pts, n, d);
    }
  }

  void delete_tree() {
    if(tree!=0) {
      delete tree;
      tree=0;
    }
  }
  List query(NumericMatrix query, const int k, const double eps=0.0) {
    return query_FR(query, k, 0.0, eps);
  }

  List query_FR(NumericMatrix query, const int k, const double radius, const double eps=0.0);

  List query_FR_ragged(NumericMatrix query, const int k, const double radius, const double eps=0.0);

  List querySelf(const int k, const double eps=0.0) {
    return queryANN_FR(data_pts, n, k, eps);
  }

  List queryWANN(const WANN& query, const int k, const double eps=0.0) {
    return queryANN_FR(query.data_pts, query.n, k, eps);
  }

  List querySelf_FR(const int k, const double radius, const double eps=0.0) {
    return queryANN_FR(data_pts, n, k, eps, radius);
  }

  List queryWANN_FR(const WANN& query, const int k, const double radius, const double eps=0.0) {
    return queryANN_FR(query.data_pts, query.n, k, eps, radius);
  }

  List queryANN_FR(const ANNpointArray query, const int nq, const int k,
                   const double eps=0.0, const double radius=0.0);

  NumericMatrix getPoints() {
    NumericMatrix points(n, d);
    for(int i = 0; i < n; i++) // Run all query points against tree
    {
      for (int j = 0; j < d; j++)
      {
        points(i,j)=data_pts[i][j];
      }
    }
    return points;
  }

};

#endif
