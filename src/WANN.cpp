#include <Rcpp.h>
#include "ANN.h"     // ANN library header
using namespace Rcpp;

RCPP_EXPOSED_CLASS(WANN)

class WANN {
  public:
  WANN(NumericMatrix data) {
    d=data.ncol();
    n=data.nrow();
    
    // now construct the points
    data_pts  = annAllocPts(n,d);
    for(int i = 0; i < n; i++) 
    {
      for(int j = 0; j < d; j++)
      {
        data_pts[i][j]=data(i,j);
      }
    }
    
    tree = new ANNkd_tree( data_pts, n, d);
  }
  
  ~WANN() {
    annDeallocPts(data_pts);
    delete tree;
  }
  
  List query(NumericMatrix query, const int k, const double eps=0.0) {
    const int nq=query.nrow();
    // ANN style point and return arrays for one point
    ANNpoint pq = annAllocPt(d);
    ANNidxArray nn_idx = new ANNidx[k];
    ANNdistArray dists = new ANNdist[k];
    
    // declare matrices for return values here
    NumericMatrix rdists(nq, k);
    IntegerMatrix ridx(nq, k);
    
    // Run all query points against tree
    for(int i = 0; i < nq; i++)
    {
      // read coords of current query point
      for(int j = 0; j < d; j++)
      {
        pq[j]=query(i,j);
      }
      
      tree->annkSearch(pq, k, nn_idx, dists, eps);
      for (int j = 0; j < k; j++)
      {
        // NB unsquare distance
        rdists(i,j) = sqrt(dists[j]);
        // put indices in returned array (nb +1 for R)
        ridx(i,j) = nn_idx[j] + 1;
      }
    }
    
    annDeallocPt(pq);
    delete [] nn_idx;
    delete [] dists;
    
    List z = List::create(Rcpp::Named("nn.idx")=ridx, Rcpp::Named("nn.dists")=rdists);
    return z ;
  }
  
  List querySelf(const int k, const double eps=0.0) {
    return queryW(this, k, eps);
  }
  
  List queryW(const WANN *query, const int k, const double eps=0.0) {
    return queryW(query->data_pts, query->n, k, eps);
  }
  
  List queryWANN(const WANN& query, const int k, const double eps=0.0) {
    return queryW(query.data_pts, query.n, k, eps);
  }
  
  List queryW(const ANNpointArray query, const int nq, const int k, const double eps=0.0) {
    
    // ANN style arrays to hold return values for one point
    ANNidxArray nn_idx = new ANNidx[k];
    ANNdistArray dists = new ANNdist[k];
    
    // declare matrices for return values here
    NumericMatrix rdists(nq, k);
    IntegerMatrix ridx(nq, k);
    
    // Run all query points against tree
    for(int i = 0; i < nq; i++) 
    {
      tree->annkSearch(query[i], k, nn_idx, dists, eps);
      for (int j = 0; j < k; j++)
      {
        // NB unsquare distance
        rdists(i,j) = sqrt(dists[j]);
        // put indices in returned array (nb +1 for R)
        ridx(i,j) = nn_idx[j] + 1;
      }
    }
    
    delete [] nn_idx;
    delete [] dists;
    
    List z = List::create(Rcpp::Named("nn.idx")=ridx, Rcpp::Named("nn.dists")=rdists);
    return z ;
  }
  
  
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
  
  private:
  ANNpointArray data_pts;
  ANNkd_tree  *tree;
  int d;
  int n;
};

RCPP_MODULE(class_WANN) {
  class_<WANN>( "WANN" )
  .constructor<NumericMatrix>()
  .method( "getPoints", &WANN::getPoints )
  .method( "query", &WANN::query )
  .method( "queryWANN", &WANN::queryWANN )
  .method( "querySelf", &WANN::querySelf )
  ;
}
