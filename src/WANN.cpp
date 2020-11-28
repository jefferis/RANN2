#include "WANN.h"

using namespace Rcpp;

WANN::WANN(NumericMatrix data, bool buildtree) : tree(0) {
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

List WANN::query_FR(NumericMatrix query, const int k, const double radius, const double eps) {
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

  const ANNdist sqRad = (ANNdist) (radius * radius);

  const bool usefr = !R_IsNA(radius) && radius > 0.0;

  // Run all query points against tree
  for(int i = 0; i < nq; i++) {

    // read coords of current query point
    for(int j = 0; j < d; j++) {
      pq[j]=query(i,j);
    }
    if(usefr) {
      tree -> annkFRSearch(
          pq, // query point
          sqRad, // squared radius
          k, // number of near neighbors to return
          nn_idx, // nearest neighbor array (modified)
          dists, // dist to near neighbors (modified)
          eps); // error bound
    } else {
      tree->annkSearch(pq, k, nn_idx, dists, eps);
    }

    // Sometimes there are fewer than k neighbors in the radius
    // In that case, set ANN's DIST_INF and NULL_ID to R's NA value
    for (int j = 0; j < k; j++) {
      if(usefr) {
        // un-square distance
        rdists(i,j) = dists[j] == ANN_DIST_INF ? NA_REAL : std::sqrt(dists[j]);
        // put indices in returned array (nb +1 for R)
        ridx(i,j) = nn_idx[j] == ANN_NULL_IDX ? NA_INTEGER : nn_idx[j] + 1;
      } else {
        rdists(i,j) = std::sqrt(dists[j]);
        ridx(i,j) = nn_idx[j] + 1;
      }
    }
  }

  annDeallocPt(pq);
  delete [] nn_idx;
  delete [] dists;

  List z = List::create(Rcpp::Named("nn.idx")=ridx, Rcpp::Named("nn.dists")=rdists);
  return z ;
}

List WANN::query_FR_ragged(NumericMatrix query, const int k, const double radius, const double eps) {
  // build tree (in case we didn't already)
  build_tree();

  const int nq=query.nrow();
  // ANN style point and return arrays for one point
  ANNpoint pq = annAllocPt(d);
  ANNidxArray nn_idx = new ANNidx[k];
  ANNdistArray dists = new ANNdist[k];

  // declare lists for return values here
  List rdists(nq);
  List ridx(nq);

  const ANNdist sqRad = (ANNdist) (radius * radius);

  // Run all query points against tree
  for(int i = 0; i < nq; i++) {

    // read coords of current query point
    for(int j = 0; j < d; j++) {
      pq[j]=query(i,j);
    }
    int npts = 0;
    npts = tree -> annkFRSearch(
      pq, // query point
      sqRad, // squared radius
      k, // number of near neighbors to return
      nn_idx, // nearest neighbor array (modified)
      dists, // dist to near neighbors (modified)
      eps); // error bound

    // creates the vectors to contain the distances and indices
    std::vector<double> rdists_element(npts);
    std::vector<int> ridx_element(npts);
    for (int j = 0; j < npts; j++) {
      // un-square distance
      rdists_element[j] = std::sqrt(dists[j]);
      // put indices in returned array (nb +1 for R)
      ridx_element[j] = nn_idx[j] + 1;
    }
    // adds the distances and indices to the output ragged lists
    rdists[i] = rdists_element;
    ridx[i] = ridx_element;
  }

  annDeallocPt(pq);
  delete [] nn_idx;
  delete [] dists;

  List z = List::create(Rcpp::Named("nn.idx")=ridx, Rcpp::Named("nn.dists")=rdists);
  return z ;
}

List WANN::queryANN_FR(const ANNpointArray query, const int nq, const int k,
                       const double eps, const double radius) {

  // build tree (in case we didn't already)
  build_tree();

  // ANN style arrays to hold return values for one point
  ANNidxArray nn_idx = new ANNidx[k];
  ANNdistArray dists = new ANNdist[k];

  // declare matrices for return values here
  NumericMatrix rdists(nq, k);
  IntegerMatrix ridx(nq, k);

  const ANNdist sqRad = (ANNdist) (radius * radius);

  const bool usefr = !R_IsNA(radius) && radius > 0.0;

  // Run all query points against tree
  for(int i = 0; i < nq; i++) {

    if(usefr) {
      tree -> annkFRSearch(
          query[i], // query point
               sqRad, // squared radius
               k, // number of near neighbors to return
               nn_idx, // nearest neighbor array (modified)
               dists, // dist to near neighbors (modified)
               eps); // error bound
    } else {
      tree -> annkSearch(query[i], k, nn_idx, dists, eps);
    }

    // Sometimes there are fewer than k neighbors in the radius
    // In that case, set ANN's DIST_INF and NULL_ID to R's NA value
    for (int j = 0; j < k; j++) {
      // un-square distance
      if(usefr) {
        rdists(i,j) = dists[j] == ANN_DIST_INF ? NA_REAL : std::sqrt(dists[j]);
        // put indices in returned array (nb +1 for R)
        ridx(i,j) = nn_idx[j] == ANN_NULL_IDX ? NA_INTEGER : nn_idx[j] + 1;
      } else {
        // NB un-square distance
        rdists(i,j) = std::sqrt(dists[j]);
        // put indices in returned array (nb +1 for R)
        ridx(i,j) = nn_idx[j] + 1;
      }
    }
  }

  delete [] nn_idx;
  delete [] dists;

  List z = List::create(Rcpp::Named("nn.idx")=ridx, Rcpp::Named("nn.dists")=rdists);
  return z ;
}

RCPP_EXPOSED_CLASS(WANN)

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
    .method( "query_FR", &WANN::query_FR )
    .method( "query_FR_ragged", &WANN::query_FR_ragged )
    .method( "queryWANN_FR", &WANN::queryWANN_FR )
    .method( "querySelf_FR", &WANN::querySelf_FR )
    ;
  }
