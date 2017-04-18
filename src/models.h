/*-----------------------------------------------------------------------------

  models.h - This is the header file for models.h
  
  Kelly Burkett; May 2007

  ---------------------------------------------------------------------------*/

#ifndef MODELS_H
#define MODELS_H

#include <vector>

#include "treeBuild.h"


using namespace std;

const double ALPHA = 2; //SNP markers only


// This function computes the prior probability of the mutation and 
// recombination rates theta and rho respectively. 
double probRates(double rate, const double low_bound, const double up_bound);


// Use an exponential/gamma prior for rho instead of a uniform
double probGamma(const double rho, const double shp=1, const double scl=0.3);


// This function computes the log of the likelihood for the intercoalescent 
// times between coalescent events min_index and max_index
double probT( vector<treeNode*> node_ptrs, int min_index, int max_index );


// The following two functions are used to compute the probability 
// of the r values. probRsided computes the value for one side 
// (left or right) for one node and probRi computes the value for both sides
// by calling probRsided
double probRsided(string side, treeNode* node_ptr, double rho,
			const vector<int>& first_markers,
			const vector<double>& locations, int ref_marker_location);
double probRi(treeNode* node_ptr, double rho,
		const vector<int>& first_markers, const vector<double>& locations,
		int ref_marker_location);
  	
	
// These two functions are used to compute the probability of the sequence
// value. Use the second version if the probability associated with the 
// sequence after the closest recombination event is not required. These
// functions call probSij and haplo/allele probs
double probSi(treeNode* node_ptr, double theta, const vector<int>& first_markers,
			 const vector<vector<double> >& allele_probs,
	     	 const vector<vector<double> >& hap_probs);
double probSiNoHap(treeNode* node_ptr, double theta, const vector<int> first_markers);


// This uses the mutation model to compute the probability of an allele given
// the children and parental alleles at that node
double probSij(string side, treeNode* node_ptr, int j, double theta,
			   const vector<int>& first_markers, const vector<vector<double> >& allele_probs,
	 		   const vector<vector<double> >& hap_probs);



// These two functions are used to compute the probabilities of the sequence
// after the closest recombination event by assuming a first-order Markov
// model for the haplotypes.
void alleleProbs(vector<string> seqs, vector<vector<double> >& allele_probs);
void haploProbs(vector<string> seqs, vector<vector<double> >& haplo_probs);

#endif
