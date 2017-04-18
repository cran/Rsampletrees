/*----------------------------------------------------------------------------

  P1P2.cpp - This is the header file for the proposal functions updating the
	         mutation and recombination rates. 

  Kelly Burkett; July 2009

  --------------------------------------------------------------------------*/

#ifndef P1P2_H
#define P1P2_H

#include <vector>

#include "nodefunctions.h"


// Proposal Distributions for either theta or rho:

// Uniform
double proposeRateUnif( bool update, double& new_rate, const double old_rate,
						const double rangefactor);

// Altered Uniform
double proposeRateUnif_2( bool update, double& new_rate, const double old_rate,
						  const double lowerbound, const double upperbound,
						  const double rangefactor );

// Normal
double proposeRateNormal( bool update, double& new_rate, const double mean,
						  const double sd );

//Gamma
double proposeRateGamma( bool update, double& new_rate, const double shape,
					     const double scale );



// P1 - This proposes a new value for theta and determines whether it is 
// accepted or rejected. A value of 0 is returned if it is rejected and 1 if it
// is accepted
int updateP1( double& theta, vector<treeNode*> node_ptrs,
			  vector<int> first_markers,
		      const double min_theta, const double max_theta );


// P2 - This proposes a new value for rho and determines whether it is 
// accepted or rejected. A value of 0 is returned if it is rejected and 1 if it
// is accepted
int updateP2( double& rho, vector<treeNode*> node_ptrs,
			  vector<int> first_markers,
			  vector<double> locations, int ref_marker_location,
			  const double min_theta, const double max_theta );

int updateP1_gamma( double& theta, vector<treeNode*> node_ptrs, vector<int> first_markers,
					 const double shape, const double scale );

int updateP2_gamma( double& rho, vector<treeNode*> node_ptrs,
				    vector<int> first_markers,
					vector<double> locations, int ref_marker_location,
					const double shape=1, const double scale=0.3 );


#endif
