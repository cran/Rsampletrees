/*----------------------------------------------------------------------------

  rng.h - this is the header file for rng.cpp

  Kelly Burkett; July 2009

  Modification Notes:

  --------------------------------------------------------------------------*/

#include <vector>
#include <gsl/gsl_rng.h>

#ifndef RNG_H
#define RNG_H

using namespace std;

// This is needed since this variable will be defined elsewhere
extern gsl_rng *rng;


unsigned long int random_seed();
		
// This function generates a U(0,1) random variable 
double generateStdUnif();


// This function generates a Uniform (lbound, ubound) random variable
double generateUnif( double lbound, double ubound );


// This function generates a discrte Uniform (lbound,ubound)  random variable
int generateDiscreteUnif( int lbound, int ubound );


// This is a generic function used to generate an exponential variable
// with a given parameter lambda ( f(x) = lambda * exp( -lambda*x ) )
double generateExp( double lambda );

// This function generates a gamma random variable with shape and
// scale parameter
double generateGamma(const double shp, const double scl );

// This generates a N(0,sd^2) RV
double generateNormal( double mean, double sd );

// This is a generic function to sample a value for a variable based on the
// probabilities in the vector probs. It returns the index of the element of
// probs chosen.
int samplePMF( vector<double> probs);


#endif



