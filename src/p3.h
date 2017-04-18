/*----------------------------------------------------------------------------

  p3.h - This is the header file for the local update functions.

  Kelly Burkett; July 2009

  Modification Notes:

  --------------------------------------------------------------------------*/

#ifndef P3_H
#define P3_H

#include <vector>

#include "nodefunctions.h"

using namespace std;

// This function is the main function to complete the local updates. Each 
// internal node is visited in turn, and at each node new values for t, r_c1,
// r_c2 and s are sampled. The MH ratio is computed and the sampled values
// are either accepted or rejected. If rejected, the old values are returned.
int updateP3(vector<treeNode*>& node_ptrs, vector<int> first_markers,
			  int ref_marker_location, vector<double> locations, double rho,
	 		  double theta, vector<vector<double> > allele_probs,
   			  vector<vector<double> > hap_probs, int i);

// This function computes the term proportional to the likelihood and is used
// to determine the acceptance probability.
double probP3( treeNode* node_ptr,
			   double theta, double rho, vector<int>& first_markers,
			   int ref_marker_location, vector<double>& locations,
			   vector<vector<double> >& allele_probs,
		       vector<vector<double> >& hap_probs,
  			   vector<treeNode*> node_ptrs );


int updateP3_allnodes( vector<treeNode*>& node_ptrs, vector<int> first_markers,
			   int ref_marker_location,
			   vector<double> locations, double rho, double theta,
			   vector<vector<double> > allele_probs,
			   vector<vector<double> > hap_probs, int it );

double probP3_all( vector<treeNode*> node_ptrs, double theta, double rho,
		vector<int>& first_markers, int ref_marker_location,
		vector<double>& locations, vector<vector<double> >& allele_probs,
		vector<vector<double> >& hap_probs );

double updateS_child(bool sample, treeNode* node_ptr,
				  const vector<int>& first_markers, double theta,
	  			  const vector<vector<double> >& allele_probs,
				  const vector<vector<double> >& hap_probs);

void getSijPMF_child(string side, treeNode* node_ptr, vector<double>& probs,
			   const vector<int>& first_markers, int locus, double theta,
	  		   const vector<vector<double> >& allele_probs,
   			   const vector<vector<double> >& hap_probs);
#endif
