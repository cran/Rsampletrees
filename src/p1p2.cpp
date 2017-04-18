/*----------------------------------------------------------------------------

  P1P2.cpp - This file contains the proposal functions for the update to the
	         mutation and recombination rates. 

  Kelly Burkett; July 2009

  Modification Notes:
  * added values[2] = exp(probRatio); for both functions since the error
    checking done assumed that log's had not been taken (Nov 2009)

  --------------------------------------------------------------------------*/

#include <vector>
#include <cmath>
#include <iostream> //For cerr/cout
#include <gsl/gsl_randist.h>
#include <Rcpp.h>

#include "nodefunctions.h"
#include "proposal.h"
#include "models.h"
#include "rng.h"
#include "p1p2.h"
#include "treeBuild.h"



// The following three functions are proposal distributions for theta and/or
// rho.

// Uniform proposal: proposed values are from U(rate/factor, rate*factor)
double proposeRateUnif( bool update, double& new_rate, const double old_rate,
						const double rangefactor )
{

	double lbound = old_rate / rangefactor;
	double ubound = rangefactor * old_rate;

	if (update==1){
		new_rate = generateUnif(lbound, ubound);
	}
	return( 1/(ubound-lbound) );
}


// Altered Uniform proposal: Don't propose values outside of the uniform prior
double proposeRateUnif_2( bool update, double& new_rate, const double old_rate,
						  const double lowerbound, const double upperbound,
						  const double rangefactor )
{

	double lbound = max(lowerbound, old_rate / rangefactor);
	double ubound = min(upperbound, rangefactor * old_rate);
	
	if (update==1){
		new_rate = generateUnif(lbound, ubound);
	}
	return( 1/(ubound-lbound) );
}


// Normal proposal: proposed values are from N(mean,sd)
double proposeRateNormal( bool update, double& new_rate, const double mean,
						  const double sd )
{

	if (update==1){
		new_rate = generateNormal( mean, sd );
	}
	return( gsl_ran_gaussian_pdf( (new_rate-mean), sd) );

}


// Gamma proposal: proposed values are from Gamma(shape, scale)
double proposeRateGamma( bool update, double& new_rate, const double shape,
							 const double scale )
{
	if (update==1){
		new_rate = generateGamma( shape, scale );
	}

	return( gsl_ran_gamma_pdf( new_rate, shape, scale) );

}



// This function completes the P1 update. It proposes a new value for theta
// from a uniform distribution. It then decides whether this value is accepted
// or rejected by computing the MH ratio and comparing that value to a uniform
// variate.
int updateP1( double& theta, vector<treeNode*> node_ptrs, vector<int> first_markers,
			  const double min_theta, const double max_theta )
{
	double new_theta = 0;
	double qA = 1, qAnew=1, MHratio = 1, probRatio = log( double(1) );
	double myunif = generateStdUnif();
	int accept=0;

	// Propose a new value and find the acceptance probability
	qAnew = proposeRateUnif_2( 1, new_theta, theta, min_theta, max_theta, 2 );
	qA = proposeRateUnif_2( 0, theta, new_theta, min_theta, max_theta, 2 );

	for (unsigned int i=0; i<node_ptrs.size()-1; i++){
		probRatio = probRatio +
				    log(probSiNoHap( node_ptrs[i], new_theta, first_markers ))-
				    log(probSiNoHap( node_ptrs[i], theta, first_markers ));
	}

	// Check that there are no infinite or 0 values as this will cause an
	// undefined MH ratio
	bool errorFlag=0;

	errorFlag = checkValue( qAnew, 1 );
	if (errorFlag==1) {
		/*cout*/ Rcpp::Rcout<< "WARNING: Probability of proposed theta is infinity" << endl;
	}

	errorFlag = checkValue( qA, 0 );
	if (errorFlag==1) {
		Rcpp::Rcout << "ERROR: Probability of current theta is 0" << endl; // cerr << "ERROR: Probability of current theta is 0" << endl;
		throw Rcpp::exception("ERROR: Probability of current theta is 0"); //exit(1);
	}

	errorFlag = checkValue( exp(probRatio), 1 );
	if (errorFlag==1) {
			/*cout*/ Rcpp::Rcout << "WARNING: Probability of data given proposed theta"
				 << "is infinity. Watch for getting stuck on this theta" << endl;
	}



	
	// Choose to accept or reject. Accept new value if the uniform number
	// generated is less than acceptance probability.
	MHratio = min( exp(probRatio+log(qA)-log(qAnew)), double( 1 ) );
	

	if ( myunif < MHratio ){
		theta = new_theta;
		accept = 1;
	}

	return( accept );
}


// This function completes the P2 update. It proposes a new value for rho
// from a uniform distribution. It then decides whether this value is accepted
// or rejected by computing the MH ratio and comparing that value to a uniform
// variate.
int updateP2( double& rho, vector<treeNode*> node_ptrs, vector<int> first_markers,
					 vector<double> locations, int ref_marker_location,
					 const double min_rho, const double max_rho )
{
	double new_rho = 0;
	double qA = 1, qAnew=1, MHratio = 1, probRatio = log(double(1));
	double myunif = generateStdUnif();
	int accept=0;


	// Propose a new value and find the probability
	qAnew = proposeRateUnif_2( 1, new_rho, rho, min_rho, max_rho, 2 );
    qA = proposeRateUnif_2( 0, rho, new_rho, min_rho, max_rho, 2);

	for (unsigned int i=0; i<node_ptrs.size()-1; i++){
		probRatio = probRatio +
					log(probRi( node_ptrs[i], new_rho, first_markers, locations,
						       ref_marker_location )) -
					log(probRi( node_ptrs[i], rho, first_markers, locations,
					           ref_marker_location ));
	}

	// Check that there are no infinite or 0 values as this will cause an
	// undefined MH ratio
	bool errorFlag=0;

	errorFlag = checkValue( qAnew, 1 );
	if (errorFlag==1) {
		/*cout*/ Rcpp::Rcout << "WARNING: Probability of proposed rho is infinity" << endl;
	}

	errorFlag = checkValue( qA, 0 );
	if (errorFlag==1) {
		Rcpp::Rcout << "ERROR: Probability of current rho is 0" << endl; // cerr << "ERROR: Probability of current rho is 0" << endl;
		throw Rcpp::exception("ERROR: Probability of current rho is 0"); //exit(1);
	}

	errorFlag = checkValue( exp(probRatio), 1 );
	if (errorFlag==1) {
			/*cout*/ Rcpp::Rcout << "WARNING: Probability of data given proposed rho"
				 << "is infinity. Watch for getting stuck on this rho" << endl;
	}

		
	// Choose to accept or reject. Accept new value if the uniform number
	// generated is less than acceptance probability.
	MHratio = min( exp(probRatio + log(qA)-log(qAnew)), double( 1 ) );

	if ( myunif < MHratio ){
		rho = new_rho;
		accept = 1;
	}

	return( accept );

}



// Using a gamma/exponential prior for theta
int updateP1_gamma( double& theta, vector<treeNode*> node_ptrs, vector<int> first_markers,
					 const double shape, const double scale )
{
	double new_theta=0;
	double qA = 1, qAnew=1, MHratio = 1, probRatio = log(double(1));
	double probA=0, probAnew=0;
	double myunif = generateStdUnif();
	int accept=0, n=0;

	n=node_ptrs.size()-1;

	//qAnew = proposeRateNormal( 1, new_theta, theta, shape*scale*scale/10 );
	//qA = proposeRateNormal( 0, theta, new_theta, shape*scale*scale/10 );
	qAnew = proposeRateUnif( 1, new_theta, theta, 2 );
	qA = proposeRateUnif( 0, theta, new_theta, 2 );


	// It is possible to generate rho values that are less than 0. If so,
	// those values are not accepted
	if (new_theta>0){

		// Determine the acceptance probability
		probA = probGamma(theta, shape, scale);
		probAnew = probGamma(new_theta, shape, scale);


		for (int i=0; i<n; i++){
			probRatio = probRatio +
						log(probSiNoHap( node_ptrs[i], new_theta, first_markers ))-
						log(probSiNoHap( node_ptrs[i], theta, first_markers ));
		}


		// Do some error checking
		bool errorFlag=0;

		errorFlag = checkValue( qAnew, 2 );
		if (errorFlag==1) {
			Rcpp::Rcout << "ERROR: Probability of proposed theta is infinity or 0" << endl; // cerr << "ERROR: Probability of proposed theta is infinity or 0" << endl;
			throw Rcpp::exception("ERROR: Probability of proposed theta is infinity or 0"); //exit(1);
		}

		errorFlag = checkValue( qA, 0 );
		if (errorFlag==1) {
			Rcpp::Rcout << "ERROR: Probability of proposed theta is infinity or 0" << endl; // cerr << "ERROR: Probability of current theta is 0" << endl;
			throw Rcpp::exception("ERROR: Probability of proposed theta is infinity or 0"); //exit(1);
		}

		errorFlag = checkValue( exp(probRatio), 1 );
		if (errorFlag==1) {
				/*cout*/ Rcpp::Rcout << "WARNING: Probability of data given proposed theta"
					 << " is infinity. Watch for getting stuck on this theta" << endl;
		}

		MHratio = min( exp(probRatio)*probAnew/probA*qA/qAnew, double( 1 ) );

	} else {
		accept = -1;
		MHratio = 0;
	}


	// Choose to accept or reject. Accept new value if the uniform number
	// generated is less than acceptance probability.
	if ( myunif < MHratio ){
		theta = new_theta;
		accept = 1;

	}

	return( accept );
}


// Using a gamma/exponential prior for rho
int updateP2_gamma( double& rho, vector<treeNode*> node_ptrs, vector<int> first_markers,
					 vector<double> locations, int ref_marker_location,
					 const double shape, const double scale )
{
	double new_rho=0;
	double qA = 1, qAnew=1, MHratio = 1, probRatio = log(double(1));
	double probA=0, probAnew=0;
	double myunif = generateStdUnif();
	int accept=0;
	unsigned int i=0;

	// Propose a new value from a normal distribution
	//qAnew = proposeRateNormal( 1, new_rho, rho, shape*scale*scale/130 );
	//qA = proposeRateNormal( 0, rho, new_rho, shape*scale*scale/130 );
	qAnew = proposeRateUnif( 1, new_rho, rho, 2 );
	qA    = proposeRateUnif( 0, rho, new_rho, 2 );

	// It is possible to generate rho values that are less than 0. If so,
	// those values are not accepted
	if (new_rho>0)
	{
		// Determine the acceptance probability
		probA = probGamma(rho, shape, scale);
		probAnew = probGamma(new_rho, shape, scale);

		for (i=0; i<node_ptrs.size()-1; i++){
			probRatio = probRatio +
						log(probRi( node_ptrs[i], new_rho, first_markers, locations,
						       ref_marker_location )) -
						log(probRi( node_ptrs[i], rho, first_markers, locations,
					           ref_marker_location ));
		}

		bool errorFlag=0;

		errorFlag = checkValue( qAnew, 2 );
		if (errorFlag==1) {
			Rcpp::Rcout << "ERROR: Probability of proposed rho is infinity or 0" << endl; // cerr << "ERROR: Probability of proposed rho is infinity or 0" << endl;
			throw Rcpp::exception("ERROR: Probability of proposed rho is infinity or 0"); //exit(1);
		}

		errorFlag = checkValue( qA, 0 );
		if (errorFlag==1) {
			Rcpp::Rcout << "ERROR: Probability of current rho is 0" << endl; // cerr << "ERROR: Probability of current rho is 0" << endl;
			throw Rcpp::exception("ERROR: Probability of current rho is 0"); //exit(1);
		}

		errorFlag = checkValue( exp(probRatio), 1 );
		if (errorFlag==1) {
				/*cout*/ Rcpp::Rcout << "WARNING: Probability of data given proposed rho"
					 << " is infinity. Watch for getting stuck on this rho" << endl;
		}

		MHratio = min( exp(probRatio)*probAnew/probA*qA/qAnew, double( 1 ) );

	} else {
		accept = -1;
		MHratio = 0;
	}


	// Choose to accept or reject. Accept new value if the uniform number
	// generated is less than acceptance probability.
	if ( myunif < MHratio ){
		rho = new_rho;
		accept = 1;
		
	}

	return( accept );
}
