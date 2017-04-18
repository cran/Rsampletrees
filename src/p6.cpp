/*----------------------------------------------------------------------------

  P6.cpp - This file contains the proposal functions for the update that 
	       swaps coalescent times for two nodes

  Kelly Burkett; July 13, 2009

  Modification Notes:

  --------------------------------------------------------------------------*/

#include <vector>
#include <iostream>
#include <Rcpp.h>


#include "nodefunctions.h"
#include "proposal.h"
#include "models.h"
#include "rng.h"
#include "p6.h"


// This function completes one step in the P6 chain. That is, it randomly
// samples two nodes from the set of compatible nodes and swaps their node
// times. The MH acceptance probability it determined and from that, whether
// to accept or reject the change is decided
int updateP6( vector<treeNode*>& node_ptrs, double theta, double rho,
			  vector<int>& first_markers, int ref_marker,
			  vector<double>& locations )
{
	int accept = 1;
	double qA = 0, qAnew = 0, MHratio = 0, probA = 0, probAnew = 0;
	double myunif = generateStdUnif();

	vector<int> row( 2, 0 );
	vector<vector<int> > compatible_mat( 0, row );
	vector<int> apair( 2 );

	// Find the set of pairs that can have their times swapped
	int numpairs = findCompatiblePairs( compatible_mat, node_ptrs );
	
	if (numpairs == 0){
		// This update can't be done so return
		return( -1 );
	}

	qAnew = 1/double(numpairs);

	
	// Select a pair from the set of pairs that can be swapped
	int pairID = generateDiscreteUnif(0, numpairs-1);
	apair = compatible_mat[pairID];

	
	// Make the change and compute the MH ratio
	probA = probP6( apair, node_ptrs, theta, rho, first_markers,
					ref_marker, locations );
	reorderEvents( apair, node_ptrs );
	qA = 1/double(findCompatiblePairs( compatible_mat, node_ptrs ));

	probAnew = probP6( apair, node_ptrs, theta, rho, first_markers,
					   ref_marker, locations );


	// Check that there are no infinite or 0 values as this will cause an
	// undefined MH ratio
	bool errorFlag=0;

	errorFlag = checkValue( qAnew, 1 );
	if (errorFlag==1) {
		/*cout*/ Rcpp::Rcout << "WARNING: Probability of proposed change in time swap update is infinity" << endl;
	}

	errorFlag = checkValue( qA, 0 );
	if (errorFlag==1) {
		/*cerr*/ Rcpp::Rcout << "ERROR: Probability of proposing current values in time swap update is 0" << endl;
		throw Rcpp::exception("ERROR: Probability of proposing current values in time swap update is 0"); //exit(1);
	}

	errorFlag = checkValue( probA, 0 );
	if (errorFlag==1) {
		/*cerr*/ Rcpp::Rcout << "ERROR: Probability of previous data is 0 in time swap update" << endl;
		throw Rcpp::exception("ERROR: Probability of previous data is 0 in time swap update"); //exit(1);
	}

	errorFlag = checkValue( probAnew, 1 );
	if (errorFlag==1) {
		/*cout*/ Rcpp::Rcout << "WARNING: Probability of data given proposed changes in time swap update"
			 << "is infinity. Watch for getting stuck" << endl;
	}

	
	// Determine the acceptance ratio and if not accepted, return to old value
	MHratio = min( probAnew/probA * qA/qAnew, double(1) );
	
	if ( myunif > MHratio ){
		reorderEvents( apair, node_ptrs );
		accept = 0;
	}

	return( accept );
}


// This function determines which pairs of internal nodes have times that are
// compatible for swapping. A matrix giving the pair IDs in a row is returned
int findCompatiblePairs( vector<vector<int> >& compatible_mat,
			 const vector<treeNode*>& node_ptrs )
{
	unsigned int i = 0, j = 0, num_pairs = 0;
	int n = (node_ptrs.size()+1)/2;
	double lb_i = 0, lb_j = 0, ub_i = 0, ub_j = 0;
	vector<int> myvec(2,0);

	compatible_mat.clear();
	for ( i = n; i < (node_ptrs.size()-1); i++ ){

		// Determine the range for the internal node of id i+1
		lb_i = max( node_ptrs[i]->child1->t, node_ptrs[i]->child2->t );
		ub_i = node_ptrs[i]->parent->t;

		for ( j = n; j < i; j++ ){

			// Determine the range for the internal node of id j+1
			lb_j = max( node_ptrs[j]->child1->t, node_ptrs[j]->child2->t );
			ub_j = node_ptrs[j]->parent->t;

			// Check that the nodes with ids i+1 and j+1 have times that make
			// these two nodes compatible for swapping
			if ( ( node_ptrs[i]->t > lb_j ) && ( node_ptrs[i]->t < ub_j ) &&
					( node_ptrs[j]->t > lb_i ) && ( node_ptrs[j]->t < ub_i ) ){
				myvec[0]=i;
				myvec[1]=j;
				compatible_mat.push_back(myvec);
				num_pairs++;
			} 
			
		}
		
	}

	return( num_pairs );
}


// This function is used to swap the event times at the selected nodes and
// to change the other associated variables: vector of node pointers and
// node ID's (which are in order of coalescent events)
// NOTE: Try using C++ swap for this?
void reorderEvents( vector<int> apair, vector<treeNode*>& node_ptrs )
{
	double temp_time = 0;
	int temp_id = 0;
	treeNode* node_1 = node_ptrs[apair[0]];
	treeNode* node_2 = node_ptrs[apair[1]];

	// Swap node times
	temp_time = node_2->t;
	changeT(node_2, node_1->t);
	changeT(node_1, temp_time);

	// Swap IDs
	temp_id = node_2->id;
	node_2->id = node_1->id;
	node_1->id = temp_id;

	// Swap elements of vector that points to each node
	node_ptrs[apair[0]] = node_2;
	node_ptrs[apair[1]] = node_1;

}


// This function computes the term proportional to the likelihood for
// the P6 update
double probP6( vector<int> apair, vector<treeNode*>& node_ptrs,
			   double theta, double rho, vector<int>& first_markers,
			   int ref_marker_location, vector<double>& locations )
{
	double prob = 1;
	int i = 0;

	treeNode* nodes[]={ node_ptrs[apair[0]], node_ptrs[apair[0]]->child1,
						 node_ptrs[apair[0]]->child2, node_ptrs[apair[1]],
						 node_ptrs[apair[1]]->child1, node_ptrs[apair[1]]->child2};

	for ( i = 0; i < 5; i++ ){
		prob = prob * probSiNoHap( nodes[i], theta, first_markers );
		prob = prob * probRi( nodes[i], rho, first_markers,
							 locations, ref_marker_location );
	}

	return( prob );

}

