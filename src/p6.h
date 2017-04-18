/*----------------------------------------------------------------------------

  P6.h - This is the header file for the proposal functions for the update 
         that swaps coalescent times for two nodes

  Kelly Burkett; July 13, 2009

  --------------------------------------------------------------------------*/

#ifndef P6_H
#define P6_H

#include <vector>

#include "nodefunctions.h"

using namespace std;

// P6 - This function completes the time swap proposal update. 
int updateP6( vector<treeNode*>& node_ptrs, double theta, double rho,
			  vector<int>& first_markers, int ref_marker,
	 		  vector<double>& locations);


// This function takes the pair of nodes stored in apair and swaps their
// node times t_i. It also swaps their ID's and swaps the elements in
// node_ptrs that point to these nodes
void reorderEvents(vector<int> apair, vector<treeNode*>& node_ptrs );


// This function determines the set of pairs that are eligible to be swapped.
// The results are stored in a k*2 matrix, with each row giving the indices of
// the the pair that can be swapped
int findCompatiblePairs( vector<vector<int> >& compatible_mat,
			 const vector<treeNode*>& node_ptrs);	


double probP6(vector<int> apair, vector<treeNode*>& node_ptrs,
			  double theta, double rho, vector<int>& first_markers,
			  int ref_marker_location, vector<double>& locations);

#endif
