/*----------------------------------------------------------------------------

  proposal.h - This is the header file for proposal.cpp
  
  Kelly Burkett; April 27, 2007

  --------------------------------------------------------------------------*/

#ifndef PROPOSAL_H
#define PROPOSAL_H

#include <vector>

#include "nodefunctions.h"
#include "readFiles.h" // For passing the Options struct


using namespace std;

// This is used to check that there are no 0 or infinity values when computing
// the MH ratio
bool checkValue(double value, int type);
		

// This function is used to update the t values for all proposals that make
// updates to the T. If sample=0, the value isn't updated but the proposal
// probability is returned.
double updateT( bool sample, treeNode* node_ptr );

			
// This function is used to update the r values for all proposals that make
// updates to the R. If sample=0, the value isn't updated but the proposal
// probability is returned.
double updateR( bool sample, treeNode* node_ptr, double rho,
			    const vector<int>& first_markers, const vector<double>& locations,
	  			int ref_location, int type,
  				treeNode* sib_ptr=NULL, treeNode* other_ptr=NULL );

// This function is used to update the s values for all proposals that make
// updates to the S. If sample=0, the value isn't updated but the proposal
// probability is returned.
double updateS( bool sample, treeNode* node_ptr,
	   		    const vector<int>& first_markers, double theta,
		  		const vector<vector<double> >& allele_probs,
				const vector<vector<double> >& hap_probs);

// This function is called by updateR and it determines the PMF for the update
// to the r value
void getRiPMF( string side, treeNode* node_ptr, vector<double>& probs,
			   vector<int>& validr, double rho, const vector<int>& first_markers,
	  		   const vector<double>& locations, int ref, int type,
		       treeNode* other_ptr, treeNode* sib_ptr );

// This function is called by updateR and it determines the PMF for the update
// to the s value
void getSijPMF(string side, treeNode* node_ptr, vector<double>& probs,
			   const vector<int>& first_markers, int j, double theta,
	 		   const vector<vector<double> >& allele_probs,
  			   const vector<vector<double> >& hap_probs);


// This function updates the temperature parameter. For this version, the
// temperature parameter corresponds to the fixed recombination rate rho
void updateI_rho( unsigned int& I, vector<double> constants,
		vector<double> lambdas, double& rho, //stringstream& ssOutFile, //ofstream& outfile,
		vector<treeNode*> node_ptrs, vector<int> first_markers,
		vector<double> locations, int ref_marker_location);

// This function computes the log of the posterior probability of all
// variables.
double probAll_relative(vector<treeNode*> node_ptrs, const double rho,
		const double theta, const vector<int>& first_markers,
		const vector<double>& locations, int ref_marker_location,
		const vector<vector<double> >& allele_probs,
		const vector<vector<double> >& hap_probs,
		const double shp=1, const double scl=1);

#endif
