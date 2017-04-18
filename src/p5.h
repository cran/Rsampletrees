/*----------------------------------------------------------------------------

  p5.h - This is the header file for the functions used to perform the
         minor topology rearrangement 
  
  Kelly Burkett; July 2009

  --------------------------------------------------------------------------*/

#ifndef P5_H
#define P5_H

#include <vector>

#include "nodefunctions.h"

using namespace std;

// This function performs the minor topology update and returns 0 if the 
// proposal was rejected and 1 if it was accepted
int updateP5(vector<treeNode*>& node_ptrs, double theta, double rho,
			  vector<int> first_markers, int ref_marker, vector<double> locations,
		   	  vector<vector<double> > hap_probs, vector<vector<double> > allele_probs);

// This function determines which nodes are eligible to be moved in this
// topology rearrangement
void findValidNodesP5(vector<treeNode*>& valid_nodes, treeNode* head_node);


// This function completes the the topology change after the nodes to be
// moved have been chosen
void topologyChangeP5( treeNode* node_ptr, treeNode* aunt_ptr, treeNode* sib_ptr );


// This determines the term proportional to the likelihood of the
// augmented data and is used to compute the acceptance probability
double probP5( treeNode* node_ptr, treeNode* aunt_ptr, treeNode* sib_ptr,
			   double theta, double rho,
			   vector<int>& first_markers, int ref_marker, vector<double>& locations,
	   		   vector<vector<double> >& hap_probs,
			   vector<vector<double> >& allele_probs);

#endif
