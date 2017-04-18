/*----------------------------------------------------------------------------

  rng.cpp - this file contains functions to generate random numbers from 
            various distributions.

  NOTES:
  - these functions were written to easily change rng libraries, so most of 
    the functions are essentially just shells that call the appropriate
    function from the library currently in use
  - GSL is currently being used
  - Checked that all the included files were actually needed (July 2009)

  Kelly Burkett; July 2009

  Modification Notes:
  * Added a random_seed function so that time is not used to seed the RNG
    as it is not adequate for cluster jobs (Nov 2009)

  --------------------------------------------------------------------------*/

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <Rcpp.h>

#include "rng.h"


unsigned long int random_seed()
{
	unsigned long int seed;
	FILE *devrandom;

	if ((devrandom = fopen("/dev/random","r")) == NULL) {
		seed = time(NULL);
		/*cout*/ Rcpp::Rcout << "Using time-based seed for RNG" << endl;
	} else {
		/*int bon = */fread(&seed,sizeof(seed),1,devrandom);
		/*cout*/ Rcpp::Rcout << "Using /dev/random/ seed for RNG"  << endl;
		fclose(devrandom);
	}

	return(seed);

}



// This function generates a U(0,1) RV. 
double generateStdUnif()
{
	return( gsl_rng_uniform(rng) );
}


// This is a generic function used to generate a uniform variable
// provided that the upper and lower bound of the range is given.
double generateUnif(double lbound, double ubound)
{
	return ( gsl_ran_flat(rng,lbound,ubound) );
}


// This is a generic function used to generate a discrete uniform variable
// in [ lbound, ubound ]. gsl_rng_uniform_int generates values in [0,numbin-1]
// If value=numbin-1=ubound-lbound+1-1, then returned value is ubound. If
// value=0 then returned value is lbound. 
int generateDiscreteUnif(int lbound, int ubound)
{
	int numbin = ubound-lbound+1;
	return( gsl_rng_uniform_int(rng,numbin)+lbound );
}


// This is a generic function used to generate an exponential
// with a given parameter lambda ( f(x) = lambda * exp( -lambda*x ) ).
// gsl_ran_exponential generates from f(x) = 1/mu exp( -1/mu*x ), where mu
// is passed to the function. So need to
// pass 1/lambda since lambda is the rate not the mean.
double generateExp(double lambda)
{
	return( gsl_ran_exponential(rng, 1/lambda) );
}

double generateGamma(const double shp, const double scl )
{
	return( gsl_ran_gamma(rng, shp, scl) );
}


double generateNormal( double mean, double sd )
{
	return( gsl_ran_gaussian(rng,sd)+mean );
}

// This is a generic function to sample a value for a variable based on the
// probabilities in the vector probs. It returns the index of the element of
// probs chosen
int samplePMF( vector<double> probs )
{
	double total = 0;
	unsigned int i = 0;
	double myunif = generateStdUnif();

	// Error if probs is empty that means that there is no distribution
	// to sample from
	if (probs.size()==0){
		/*cerr*/ Rcpp::Rcout << "ERROR: The vector of probabilities to sample from "
		     << "is empty" 
		     << endl;
		return(-1); //exit(1);
	}

	// Replace probs with a cumulative distribution
	// NOTE: Check before and after values
	for ( i = 0; i < probs.size(); i++ ){
		total = probs[i] + total;
		probs[i] = total;
	}
		
	// Ensure that the CDF sums to 1. If the condition were "probs[n]==1" then
	// it would always fail because the answer is never exactly 1, it will be
	// 0.99999995 for example. So that is why the tolerance is used. There may
	// be an approximately equal function in C++ but I haven't yet found it.
	double	tol = 0.00000001;
	if (abs(1-probs[probs.size()-1])>tol){
		/*cerr*/ Rcpp::Rcout << "ERROR: PMF does not sum to 1" << endl;
		return(-1); //exit(1);
	}

	// Find which interval of the CDF of probs that myunif lies in. If it
	// can't be found there is a problem with the sampling so exit.
	i = 0;
	while ( (i<probs.size()) && (myunif > probs[i]) ){
		i++;
	}

	if (i>=probs.size()){
		/*cerr*/ Rcpp::Rcout << "ERROR: Problem in sampling the PMF" << endl;
		return(-1); //exit(1);
	}

	return( i );
}
