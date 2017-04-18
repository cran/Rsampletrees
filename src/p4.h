/*----------------------------------------------------------------------------

  p4.h - This is the header file for the functions used to perform the
         major topology rearrangement 
  
  Kelly Burkett; July 2009

  --------------------------------------------------------------------------*/

#ifndef P4_H
#define P4_H

#include <vector>

#include "nodefunctions.h"

using namespace std;

// This function performs the major topology update and returns 0 if the 
// proposal was rejected and 1 if it was accepted
int updateP4(vector<treeNode*>& node_ptrs, double theta, double rho,
			  vector<int>& first_markers, int ref_marker, vector<double>& locations,
	 		  vector<vector<double> >& hap_probs, vector<vector<double> >& allele_probs );

int updateP4_2( vector<treeNode*>& node_ptrs, double theta, double rho,
			  vector<int>& first_markers, int ref_marker,
	 		  vector<double>& locations, vector<vector<double> >& hap_probs,
	  		  vector<vector<double> >& allele_probs, int reps);


// This function determines which nodes can be moved from their current 
// location in the tree
void findValidNodesP4(vector<treeNode*>& valid_nodes, treeNode* head_node);


// Given a node has been chosen to be moved, this function determines the set
// of nodes that can accept that node as its sibling
void selectOtherNodeP4(vector<double>& probs, vector<treeNode*>& values,
				  vector<treeNode*>& node_ptrs, treeNode* chosen_node,
	 			  vector<int>& first_markers);	  


void selectOtherNodeP4_2(vector<double>& probs, vector<treeNode*>& values,
					   const vector<treeNode*>& node_ptrs, treeNode* chosen_node,
					   const vector<int>& first_markers, double theta, double rho,
					   const vector<double>& locations);

// This function implements the topology change 
void topologyChangeP4( treeNode* node_ptr, treeNode* other_ptr, treeNode* sib_ptr,
					   vector<treeNode*>& node_ptrs );


// This determines the term proportional to the log of the likelihood of the
// augmented data and is used to compute the acceptance probability
double probP4(treeNode* node_ptr, treeNode* other_ptr, treeNode* sib_ptr,
			  vector<treeNode*> node_ptrs, 
	 		  double theta, double rho, vector<int>& first_markers,
 			  int ref_marker, vector<double>& locations,
	  		  vector<vector<double> >& allele_probs, vector<vector<double> >& hap_probs,
	  		  bool outflag=0);

#endif

