
#include <Rcpp.h>
#include "ANN.h"     // ANN library header

using namespace Rcpp;

// [[Rcpp::export]]
List nn2(NumericMatrix data, NumericMatrix query, const int k, const double eps=0.0) {
	const int d=data.ncol();
	const int nd=data.nrow();
	const int nq=query.nrow();
	
	ANNkd_tree	*the_tree;	// Search structure

	ANNpointArray data_pts 	= annAllocPts(nd,d);		// Allocate data points
	ANNidxArray nn_idx 		= new ANNidx[k];		// Allocate near neigh indices
	ANNdistArray dists 		= new ANNdist[k];		// Allocate near neighbor dists
	// now construct the points
	for(int i = 0; i < nd; i++) 
	{
		for(int j = 0; j < d; j++)
		{
			data_pts[i][j]=data(i,j);
		}
	}
	the_tree = new ANNkd_tree( data_pts, nd, d);
	
	// return values here
	NumericMatrix rdists(nq, d);
	IntegerVector ridx(nq, d);
	
	//now iterate over query points
	ANNpoint pq = annAllocPt(d);
	
	for(int i = 0; i < nq; i++)	// Run all query points against tree
	{
		// read coords of current query point
		for(int j = 0; j < d; j++)
		{
			pq[j]=query(i,j);
		}
		
		the_tree->annkSearch(	// search
		pq,	// query point
		k,		// number of near neighbors
		nn_idx,		// nearest neighbors (returned)
		dists,		// distance (returned)
		eps);	// error bound
	}

	annDeallocPt(pq);
	annDeallocPts(data_pts);
	delete the_tree;
	delete [] nn_idx;
	delete [] dists;

	List z = List::create( rdists, ridx);
	return z ;
}
