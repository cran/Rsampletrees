/*-----------------------------------------------------------------------------

  models.cpp - This file contains the functions for computing the 
               probabilities of the data ( Pr(A|G) ). They are called when
               calculating the Metropolis-Hastings ratio, after a particular
               proposal has been generated from the proposal distribution

  Kelly Burkett; May 2007

  Modification Notes:
  * Aug/Sept 2009 - Flipped vector/matrix orientation of allele and
    haplo_probs for easier use when haplotypes are unknown

  ---------------------------------------------------------------------------*/

#include <cmath>
#include <bitset>
#include <iostream> // For error messages
#include <fstream>
#include <vector>
#include <gsl/gsl_randist.h>
#include <Rcpp.h>

#include "models.h"


// This function computes the prior probability of the mutation and 
// recombination rates theta and rho respectively. 
// The prior is assumed to be uniform on [low_bound, up_bound]
double probRates(double rate, const double low_bound, const double up_bound)
{
	double lbound=1, ubound=1;
	
	lbound = max(low_bound,rate/2);
	ubound = min(up_bound,2*rate);
	return( 1/(ubound-lbound) );
	
}


double probGamma(const double rho, const double shp, const double scl)
{
	if (rho>0){
		return ( gsl_ran_gamma_pdf(rho, shp, scl) );
	} else {
		return ( 0 );
	}
}



// Compute the probability of the intercoalescence times between nodes
// corresponding to min_index and max_index in the vector of node_ptrs.
// The first coalescence time has t_min_index as the lower bound. The
// last coalescence time has t_max_index as the upper bound since
// the final i is 1 less than max_index but the upper value is [i+1]
// NOTE: This is logged already.
double probT( vector<treeNode*> node_ptrs, int min_index, int max_index  )
{
	double prob=log(double(1));
	int i=0;
	double omega=0;
	int num_seq=0, total_seq=(node_ptrs.size()+1)/2;
	double expparam=0;

	
	for (i=min_index; i<=max_index; i++){
		omega = node_ptrs[i]->t - node_ptrs[i-1]->t;
		if (omega==0){
			prob=0;
		}
		num_seq = 2*total_seq-i;
		expparam = (num_seq) * (num_seq-1)/2;
		prob = prob + log(expparam*exp(-expparam*omega));
	}

	return( prob );
}


// This function implements the ZP probability for r_i.
// It is from equation A3
double probRsided(string side, treeNode* node_ptr, double rho,
				const vector<int>& first_markers, const vector<double>& locations,
			    int ref_marker_location)
{
	
	double prob=1;
	double lambda=node_ptr->b*rho/2; 
	int r=0, z=0, d1_index=0, d0_index=0, d0=0, d1=0;
	
	if ( side == "left" ){
		r = node_ptr->r_L;
		z = node_ptr->z_L;

		//Figure out the markers on either side of the recombination
		//r is to the left and r+1 is to the right.
		d1_index = first_markers[0]-r+1;
		d0_index = d1_index +1;

		//Determine the distances from the focal point to each marker
		if ( r!=1 ){
			d0=ref_marker_location-locations[d0_index];}
		if (r!=z ){
			d1=ref_marker_location-locations[d1_index];}
			
	} else {
		r = node_ptr->r_R;
		z = node_ptr->z_R;

		//Figure out the markers on either side of the recombination
		//r is to the left and r+1 is to the right
		d1_index = first_markers[1]+r-1;
		d0_index = d1_index -1;

		//Determine the distances from the focal point to each marker
		if ( r!=1 ){
			d0=locations[d0_index]-ref_marker_location;}
		if (r!=z ){
			d1=locations[d1_index]-ref_marker_location;}
		//if ( (r==1)&&(r==z) ){
		//	d0=locations[d1_index]-ref_marker_location;}
	}


	if ( r == 0 ) {
		// either all markers are to the left or right of the focal point
		// or we are at the head node. In both cases, there is no
		// probability to compute so return 1 so that the product
		// is multiplied by 1.
		prob = 1;
	} else if ( (r == 1)&&(r!=z) ){
		// Recombination between reference locus and first marker so d0=0
		prob = 1-exp( -lambda*d1);
	} else if ( r < z ){
		// Recombination between the (r-1)th and rth locus
		prob = exp(-lambda*d0)-exp(-lambda*d1);
	} else if (r == z){
		// Recombination after the last locus that we are tracking
		prob = exp(-lambda*d0);
	} else {
		// This shouldn't happen based on how we generate values so error
		prob = 0;
		//cerr << "WARNING: r > z; Problem in generating probabilities for r" << endl;
		Rcpp::Rcout << "WARNING: r > z; Problem in generating probabilities for r" << endl;
	}
	
	return(prob);
	
} 


// This function computes the probability of the recombination variables
// for a node and is used for computing pr( A | G ). It is fairly basic but
// done this way it mirrors the sequence probabilities
double probRi( treeNode* node_ptr, double rho,
			 const vector<int>& first_markers, const vector<double>& locations,
			 int ref_marker_location)
{
	double prob=1;
	prob = probRsided("left", node_ptr, rho, first_markers, locations,
				  ref_marker_location);
	prob = prob * probRsided("right", node_ptr, rho, first_markers, locations,
						 ref_marker_location);

	return(prob);
}


// This function computes the probability of the sequence data for a node
// and is used for computing pr( A | G ).
// NOTE: We could still rearrange this and probSij as a loop from z_L to z_R
//       (excluding the focal point locus if it happens to be one of the
//       markers). We still will need to convert between j and locus though
//       because of how r and z are stored relative to the focal point
double probSi( treeNode* node_ptr, double theta, const vector<int>& first_markers,
			  const vector<vector<double> >& allele_probs,
	 		  const vector<vector<double> >& hap_probs)
{
	double prob=1;
	int locus=0, z=0;

	// Get probability for sequence to the left of focal point
	z = first_markers[0]-node_ptr->z_L+1;
	for ( locus=first_markers[0]; locus>z; locus--){
		prob = prob * probSij("left", node_ptr, locus, theta, first_markers,
							  allele_probs, hap_probs);
	}

	// Get probability for sequence to the right of focal point
	z = first_markers[1]+node_ptr->z_R-1;
	for ( locus=first_markers[1]; locus<z; locus++){
		prob = prob * probSij("right", node_ptr, locus, theta, first_markers,
 							  allele_probs, hap_probs);
	} 

	return(prob);
}


// This function computes the probability of the sequence data for a node
// and is used for computing pr( A | G ). Note that we are not computing
// the probability of the haplotype portion as that will cancel for the
// particular update made. 
double probSiNoHap( treeNode* node_ptr, double theta, const vector<int> first_markers)
{
	double prob=1;
	int locus=0, r=0;
	vector<vector<double> > allele_probs(0);
	vector<vector<double> > hap_probs(0);

	// Get probability for sequence to the left of focal point
	r = first_markers[0]-node_ptr->r_L+1;
	for ( locus=first_markers[0]; locus>r; locus--){
		prob = prob * probSij("left", node_ptr, locus, theta, first_markers,
							  allele_probs, hap_probs);
	}

	// Get probability for sequence to the right of focal point
	r = first_markers[1]+node_ptr->r_R-1;
	for ( locus=first_markers[1]; locus<r; locus++){
		prob = prob * probSij("right", node_ptr, locus, theta, first_markers,
							  allele_probs, hap_probs);
	}

	return(prob);
}


// This computes the allele probability associated with a single locus. We pass
// the side and the value of j (1 to z-1). The probability comes from either
// the mutation model or the haplotype model, depending on j and z
double probSij(string side, treeNode* node_ptr, int locus, double theta,
			   const vector<int>& first_markers, const vector<vector<double> >& allele_probs,
	 		   const vector<vector<double> >& hap_probs)
{
	double probA=1;
	int prev_col=-1, prev_hcol=1;
	string prev_haplo(1,'0');
	int prev_row=0, cur_row=0;
	double lambda = node_ptr->b*theta/2;
	int r=0, type=0;
	

	// Values for the variables depends on whether we are on the left
	// or right side of the focal point
	if ( side == "left" ){
		
		r = first_markers[0]-node_ptr->r_L+1; // index of 1st marker after recomb point

		// Determine whether this locus was inherited from this node's
		// parent (type=1) or not. If not, are we at the first locus to
		// recombine in (type=2) or not (type=3).
		if ( locus > r ){
			type = 1;
		} else {
			type = 2;
			prev_col = locus+1;
			prev_hcol = locus;
			if ( ((node_ptr->parent == NULL)&&(locus != first_markers[0]))||
				((node_ptr->parent != NULL)&&(locus < r)) ) {
				type = 3;
				prev_haplo = node_ptr->seq.substr(locus,2);
			}
		}

			
	} else {

		r = first_markers[1]+node_ptr->r_R-1;

		// These are needed if this locus was not inherited from this
		// node's parent (ie it recombined in)
		if ( locus < r ){
			type = 1;
		} else {
			type = 2;
			prev_col = locus-1;
			prev_hcol = locus-1;
			if ( ((node_ptr->parent == NULL)&&(locus != first_markers[1]))||
							((node_ptr->parent != NULL)&&(locus > r)) ){
				type = 3;
				prev_haplo = node_ptr->seq.substr(locus-1,2);
			}
		}

	}

	// This converts the two-locus haplotype sequence into a binary
	// variable. This avoids the use of if statements for "00", "01" etc.
	// Note that we must specify a default for next_haplo since it is in the
	// if statement. However, this variable is only used under the conditions
	// of the if, so if the value is not changed from the default, it doesn't
	// matter
	bitset<2> p_haplo (prev_haplo);


	// This sets the variables that index which rows of the allele/haplo
	// probability matrices. Converting in one step from a char element in
	// a string to a number ended up not working as expected.
	if (node_ptr->seq[prev_col] == '0'){
		prev_row = 0;
	} else { prev_row = 1; }
	if (node_ptr->seq[locus] == '0'){
		cur_row = 0;
	} else { cur_row = 1; }


	// Determine the prob of the given locus based on how it was inherited
	// from its parent

	if ( type == 1){ // Probability computed from mutation model
		if (node_ptr->seq[locus] == node_ptr->parent->seq[locus]){
			probA = 1/ALPHA*(1-exp(-lambda))+exp(-lambda);
		} else {
			probA = 1/ALPHA*(1-exp(-lambda));
		}
	} else if ( type == 2){ //Probability computed from haplo model
		probA=allele_probs[locus][cur_row];
		//probA=allele_probs[cur_row][locus];

	} else {
		probA = hap_probs[prev_hcol][p_haplo.to_ulong()]/
				allele_probs[prev_col][prev_row];
		//probA = hap_probs[p_haplo.to_ulong()][prev_hcol]/
		//allele_probs[prev_row][prev_col];
	}

	return(probA);
	
}


// This sets up a 2-dim vector of allele probabilities for each of the
// loci. There are 2 columns (0 and 1) and seq_length rows
void alleleProbs(vector<string> seqs, vector<vector<double> >& allele_probs)
{
	unsigned int i=0, j=0;
	unsigned int nloci = seqs[0].size();
	unsigned int nseq = seqs.size();
	unsigned int nallele = 2;

	vector<double> row ( nallele, 0 );

	for (i=0; i<nloci; i++){
		allele_probs.push_back(row);
	}

	// Count up the numbers with each of the allelic types. 
	for (i=0; i<nseq; i++){ // over all sequences
		for (j=0; j<nloci; j++){ // over all loci in a sequence
			
			// I would like this to be done with a for loop but I couldn't
			// get the seqs[i][j] to convert to an integer properly. Look
			// into this.
			if (seqs[i][j] == '0'){
				allele_probs[j][0] += 1;
			}  else if (seqs[i][j] == '1'){
				allele_probs[j][1] += 1;
			} 	else {
					//cerr << "Improper allele type" << endl;
					Rcpp::Rcout << "Improper allele type" << endl;
					//exit(1);
					throw Rcpp::exception("Improper allele type");

			}
		}
	}
	
	// Now convert to proportions as described in ZP
	for (i=0;i<nloci;i++){
		for (j=0; j<nallele; j++){
			allele_probs[i][j] = (allele_probs[i][j]+1)/(nseq+2);
		}
	}
		
}



// This finds the 2-locus haplotype proportions along the sequence. The
// first column corresponds to 00, the second 01, the third 10 and the fourth
// 11. The first row corresponds to haplotype of elements 0-1 of the
// sequence vector and so on.
void haploProbs(vector<string> seqs, vector<vector<double> >& haplo_probs)
{
	unsigned int i=0, j=0;
	string haplo;
	vector<double> row(4,0);
	for (i=0; i<seqs[0].length()-1; i++){
		haplo_probs.push_back(row);
	}
	
	// Count up the numbers with each of the haplotypic types. 
	for (i=0; i<seqs.size(); i++){ // over all sequences
		for (j=0; j<(seqs[i].length()-1); j++){ // over all loci in a sequence
			haplo = seqs[i].substr(j,2);
			if (haplo == "00"){
				haplo_probs[j][0] += 1;
			}
			if (haplo == "01"){
				haplo_probs[j][1] += 1;
			}
			if (haplo == "10"){
				haplo_probs[j][2] += 1;
			}
			if (haplo == "11"){
				haplo_probs[j][3] += 1;
			}
		}
	}
	
	// Now convert to proportions as described in ZP
	for (i=0;i<haplo_probs.size();i++){
		for (j=0;j<haplo_probs[0].size();j++){
			haplo_probs[i][j] = (haplo_probs[i][j]+1)/(seqs.size()+4);
		}
	}
	
}







