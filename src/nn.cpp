
#include <Rcpp.h>
#include <ANN/ANN.h>     // ANN library header

using namespace Rcpp;

// [[Rcpp::export]]
List nn2_cpp(NumericMatrix data, NumericMatrix query, const int k, const double eps = 0.0, const double radius = NA_REAL) {

	const int d=data.ncol();
	const int nd=data.nrow();
	const int nq=query.nrow();

	ANNkd_tree	*the_tree;	// Search structure

	ANNpointArray data_pts 	= annAllocPts(nd,d);		// Allocate data points
	ANNidxArray nn_idx 		= new ANNidx[k];		// Allocate near neigh indices
	ANNdistArray dists 		= new ANNdist[k];		// Allocate near neighbor dists

	// now construct the points
	for(int i = 0; i < nd; i++) {
		for(int j = 0; j < d; j++) {
			data_pts[i][j]=data(i,j);
		}
	}

	the_tree = new ANNkd_tree( data_pts, nd, d);

	// return values here
	NumericMatrix rdists(nq, k);
	IntegerMatrix ridx(nq, k);

	//now iterate over query points
	ANNpoint pq = annAllocPt(d);


	// By default, we do a normal k nearest neighbor search
	if(R_IsNA(radius)) {

	  // Run all query points against tree
	  for(int i = 0; i < nq; i++)	{

	    // read coords of current query point
	    for(int j = 0; j < d; j++) {
	      pq[j]=query(i,j);
	    }

	    the_tree -> annkSearch(	// search
	        pq,	// query point
	        k,		// number of near neighbors
	        nn_idx,		// nearest neighbors (returned)
	        dists,		// distance (returned)
	        eps // error bound
	    );

	    for (int j = 0; j < k; j++) {
	      rdists(i,j) = sqrt(dists[j]);	// un-square the distance
	      ridx(i,j) = nn_idx[j] + 1;	// put indices in returned array (nb +1 for R)
	    }
	  }

	} else {
	  // if radius is set to some value other than R's NA, we perform
	  // a fixed radius k nearest neighbor search

	  // Run all query points against tree
	  for(int i = 0; i < nq; i++)	{

	    // read coords of current query point
	    for(int j = 0; j < d; j++) {
	      pq[j]=query(i,j);
	    }

	    const ANNdist sqRad = (ANNdist) (radius * radius);
	    the_tree -> annkFRSearch(
	        pq, // query point
	        sqRad, // squared radius
	        k, // number of near neighbors to return
	        nn_idx, // nearest neighbor array (modified)
	        dists, // dist to near neighbors (modified)
	        eps // error bound
	    );

	    // Sometimes there are fewer than k neighbors in the radius
	    // In that case, set ANN's DIST_INF and NULL_ID to R's NA value
	    for (int j = 0; j < k; j++) {
	      rdists(i,j) = dists[j] == ANN_DIST_INF ? NA_REAL : sqrt(dists[j]);	// un-square distance
	      ridx(i,j) = nn_idx[j] == ANN_NULL_IDX ? NA_INTEGER : nn_idx[j] + 1;	// put indices in returned array (nb +1 for R)
	    }
	  }
	}

	annDeallocPt(pq);
	annDeallocPts(data_pts);
	delete the_tree;
	delete [] nn_idx;
	delete [] dists;

	List z = List::create(Rcpp::Named("nn.idx")=ridx, Rcpp::Named("nn.dists")=rdists);
	return z ;
}
