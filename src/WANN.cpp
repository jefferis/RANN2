#include <Rcpp.h>
#include "ANN.h"     // ANN library header
using namespace Rcpp;

class WANN {
	public:
	WANN(NumericMatrix data) {
		d=data.ncol();
		nd=data.nrow();
		
		data_pts 	= annAllocPts(nd,d);		// Allocate data points
		
		// now construct the points
		for(int i = 0; i < nd; i++) 
		{
			for(int j = 0; j < d; j++)
			{
				data_pts[i][j]=data(i,j);
			}
		}
		the_tree = new ANNkd_tree( data_pts, nd, d);
	}
	~WANN() {
		annDeallocPts(data_pts);
		delete the_tree;
	}
	
	List query(NumericMatrix query, const int k, const double eps=0.0) {
		const int nq=query.nrow();
		ANNidxArray nn_idx 		= new ANNidx[k];		// Allocate near neigh indices
		ANNdistArray dists 		= new ANNdist[k];		// Allocate near neighbor dists
		// return values here
		NumericMatrix rdists(nq, k);
		IntegerMatrix ridx(nq, k);
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
			for (int j = 0; j < k; j++)
			{
				rdists(i,j) = sqrt(dists[j]);	// unsquare distance
				ridx(i,j) = nn_idx[j] + 1;	// put indices in returned array (nb +1 for R)
			}
		}
		
		annDeallocPt(pq);
		delete [] nn_idx;
		delete [] dists;
		
		List z = List::create(Rcpp::Named("nn.idx")=ridx, Rcpp::Named("nn.dists")=rdists);
		return z ;
	}
	
	private:
		ANNpointArray data_pts;
		ANNkd_tree	*the_tree;
		int d;
		int nd;
};

RCPP_MODULE(class_WANN) {
    class_<WANN>( "WANN" )
    .constructor<NumericMatrix>()
    .method( "query", &WANN::query )
    ;
}
