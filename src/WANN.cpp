#include <Rcpp.h>
#include <ANN/ANN.h>     // ANN library header
using namespace Rcpp;

RCPP_EXPOSED_CLASS(WANN)

class WANN {
  public:
  WANN(NumericMatrix data, bool buildtree=true) : tree(0) {
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

    if(buildtree) build_tree();
  }

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

  List query(NumericMatrix query, const int k, const double eps=0.0, const double radius = NA_REAL) {
    // build tree (in case we didn't already)
    build_tree();

    const int nq=query.nrow();
    // ANN style point and return arrays for one point
    ANNpoint pq = annAllocPt(d);
    ANNidxArray nn_idx = new ANNidx[k];
    ANNdistArray dists = new ANNdist[k];

    // declare matrices for return values here
    NumericMatrix rdists(nq, k);
    IntegerMatrix ridx(nq, k);

    // By default, we do a normal k nearest neighbor search
    if(R_IsNA(radius)) {

      // Run all query points against tree
      for(int i = 0; i < nq; i++) {

        // read coords of current query point
        for(int j = 0; j < d; j++) {
          pq[j]=query(i,j);
        }

        tree->annkSearch(pq, k, nn_idx, dists, eps);

        for (int j = 0; j < k; j++) {
          // NB unsquare distance
          rdists(i,j) = sqrt(dists[j]);
          // put indices in returned array (nb +1 for R)
          ridx(i,j) = nn_idx[j] + 1;
        }
      }
    } else {
      // if radius is set to some value other than R's NA, we perform
      // a fixed radius k nearest neighbor search

      const ANNdist sqRad = (ANNdist) (radius * radius);

      // Run all query points against tree
      for(int i = 0; i < nq; i++) {

        tree -> annkFRSearch(
            pq, // query point
            sqRad, // squared radius
            k, // number of near neighbors to return
            nn_idx, // nearest neighbor array (modified)
            dists, // dist to near neighbors (modified)
            eps); // error bound

        // Sometimes there are fewer than k neighbors in the radius
        // In that case, set ANN's DIST_INF and NULL_ID to R's NA value
        for (int j = 0; j < k; j++) {
          rdists(i,j) = dists[j] == ANN_DIST_INF ? NA_REAL : std::sqrt(dists[j]);	// un-square distance
          ridx(i,j) = nn_idx[j] == ANN_NULL_IDX ? NA_INTEGER : nn_idx[j] + 1;	// put indices in returned array (nb +1 for R)
        }
      }
    }

    annDeallocPt(pq);
    delete [] nn_idx;
    delete [] dists;

    List z = List::create(Rcpp::Named("nn.idx")=ridx, Rcpp::Named("nn.dists")=rdists);
    return z ;
  }

  List querySelf(const int k, const double eps=0.0, const double radius = NA_REAL) {
    return queryANN(data_pts, n, k, eps, radius);
  }

  List queryWANN(const WANN& query, const int k, const double eps=0.0, const double radius = NA_REAL) {
    return queryANN(query.data_pts, query.n, k, eps, radius);
  }

  List queryANN(const ANNpointArray query, const int nq, const int k, const double eps=0.0, const double radius = NA_REAL) {

    // build tree (in case we didn't already)
    build_tree();

    // ANN style arrays to hold return values for one point
    ANNidxArray nn_idx = new ANNidx[k];
    ANNdistArray dists = new ANNdist[k];

    // declare matrices for return values here
    NumericMatrix rdists(nq, k);
    IntegerMatrix ridx(nq, k);

    // By default, we do a normal k nearest neighbor search
    if(R_IsNA(radius)) {

      // Run all query points against tree
      for(int i = 0; i < nq; i++) {

        tree -> annkSearch(query[i], k, nn_idx, dists, eps);

        for (int j = 0; j < k; j++) {
          // NB un-square distance
          rdists(i,j) = std::sqrt(dists[j]);
          // put indices in returned array (nb +1 for R)
          ridx(i,j) = nn_idx[j] + 1;
        }
      }

    } else {

      // if radius is set to some value other than R's NA, we perform
      // a fixed radius k nearest neighbor search

      const ANNdist sqRad = (ANNdist) (radius * radius);

      // Run all query points against tree
      for(int i = 0; i < nq; i++) {

        tree -> annkFRSearch(
            query[i], // query point
                 sqRad, // squared radius
                 k, // number of near neighbors to return
                 nn_idx, // nearest neighbor array (modified)
                 dists, // dist to near neighbors (modified)
                 eps); // error bound

        // Sometimes there are fewer than k neighbors in the radius
        // In that case, set ANN's DIST_INF and NULL_ID to R's NA value
        for (int j = 0; j < k; j++) {
          rdists(i,j) = dists[j] == ANN_DIST_INF ? NA_REAL : std::sqrt(dists[j]);	// un-square distance
          ridx(i,j) = nn_idx[j] == ANN_NULL_IDX ? NA_INTEGER : nn_idx[j] + 1;	// put indices in returned array (nb +1 for R)
        }
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
  .constructor<NumericMatrix,bool>()
  .method( "getPoints", &WANN::getPoints )
  .method( "build_tree", &WANN::build_tree )
  .method( "delete_tree", &WANN::delete_tree )
  .method( "query", &WANN::query )
  .method( "queryWANN", &WANN::queryWANN )
  .method( "querySelf", &WANN::querySelf )
  ;
}
