/*----------------------------------------------------------------------------

  proposal.cpp - This file contains functions that are used in making a subset
		 of the proposals 1-6. Since multiple updates access these
		 functions, they are stored in a separate file. 

  Kelly Burkett; April 2007

  Modification Notes:
  * rearranged project so that the functions for P1 - P6 were in their own
    files. This one now contains those functions called by P1 - P6 (July 2009)

  --------------------------------------------------------------------------*/

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <Rcpp.h>

#include "nodefunctions.h"
#include "models.h"
#include "rng.h"
#include "misc.h"
#include "proposal.h"


// This function is used to check if the value passed is either infinity
// or 0. If so, this may cause problems so return a flag indicating
// that this has occurred
bool checkValue( double value, int type)
{
	bool flag=0;
			
	if ( (type==0)||(type==2) ){
		if ( value == 0 ){
				/*cout*/ Rcpp::Rcout << "WARNING: 0 value" << endl;
				flag=1;
		}
	}
	if ( (type==1)||(type==2) ){
		if ( value == numeric_limits<double>::infinity() ){
				/*cout*/ Rcpp::Rcout << "WARNING: infinity value" << endl;
				flag=1;
		}
	}
	
	return(flag);
}



// This function is used to update the time for the local updates and the
// major topology rearrangement. The time is generated from a uniform
// distribution if the node is not the head node, and from an exponential(1)
// distribution if it is the head node.
double updateT( bool sample, treeNode* node_ptr )
{
	
	double tcmax = max( node_ptr->child1->t, node_ptr->child2->t );

	double t=0, prob=1;

	// Sample a new t. This is done differently depending on whether this is
	// the head node or not
	if ( node_ptr->parent == NULL ) {
		
		//This is the head node, waiting time exp(1)
		if ( sample == 1){
			t = generateExp( 1 );
			if (t==0){
				/*cout*/ Rcpp::Rcout << "WARNING: A value of 0 was generated, from exp()" << endl;
			}
			changeT(node_ptr, t+tcmax);
		} else {
			t = node_ptr->t - tcmax;
		}
		prob = exp(-t);
		if (prob==1){
			/*cout*/ Rcpp::Rcout << "prob of time is 1 in updateT" << endl;
		}
	
	} else {

		// t has a U(tcmax,ta) distribution
		double ta = node_ptr->parent->t;
		
		if ( sample == 1){
			t = generateUnif(tcmax, ta);
			if (t==0){
				/*cout*/ Rcpp::Rcout << "WARNING: A value of 0 was generated, from unif" << endl;
				/*cout*/ Rcpp::Rcout << "tcmax=" << tcmax << "; ta=" << ta << endl;
			}
			changeT(node_ptr, t);
		}
		prob = 1/(ta-tcmax);

	}
	
	return( prob );

}


// This  function is used for the proposal portion associated with new
// recombination variables. The PMF is determined for all possible values of
// r given z and other information in the tree. Then either a new value
// for that locus is sampled based on the PMF and the probability is returned
// or the probability of the actual value, with respect to the PMF is returned.
double updateR(bool sample, treeNode* node_ptr, double rho,
			   const vector<int>& first_markers, const vector<double>& locations,
			   int ref_location, int type, treeNode* other_ptr,
			   treeNode* sib_ptr)
{
	double prob=1;
	vector<double> probs( 0 );
	vector<int> rvalues( 0 );
	int r_index=0;
	vector<int>::iterator it;


	// Deal with left side
	getRiPMF("left", node_ptr, probs, rvalues, rho, first_markers, locations,
			 ref_location, type, other_ptr, sib_ptr);
	
	if (sample == 1){
		r_index = samplePMF( probs );
		changeR( "left", node_ptr,  rvalues[r_index] );
		prob = prob * probs[r_index];
	} else {
		it = find (rvalues.begin(), rvalues.end(), node_ptr->r_L);
		if ( it == rvalues.end() ){
			/*cerr*/ Rcpp::Rcout << "ERROR: rvalue is not valid" << endl;
			/*cerr*/ Rcpp::Rcout << "To find: " << node_ptr->r_L << " " << probs[int(it-rvalues.begin())] << endl;
			throw Rcpp::exception("ERROR: rvalue is not valid To find:"); //exit (1);
		} else {
			prob = prob* probs[int(it-rvalues.begin())];
		}
	}
   

	// Deal with the right side
	getRiPMF("right", node_ptr, probs, rvalues, rho, first_markers, locations,
			 ref_location, type, other_ptr, sib_ptr);

	if (sample == 1){
		r_index = samplePMF( probs );
		changeR( "right", node_ptr, rvalues[r_index] );
		prob = prob * probs[r_index];
	} else {
		it = find (rvalues.begin(), rvalues.end(), node_ptr->r_R );
		if ( it == rvalues.end() ){
			/*cerr*/ Rcpp::Rcout << "ERROR: rvalue is not valid" << endl;
			/*cerr*/ Rcpp::Rcout << "To find: " << node_ptr->r_L << " "
				 <<	probs[int(it-rvalues.begin())] << endl;
			throw Rcpp::exception("ERROR: rvalue is not valid"); //exit (1);
		} else {
			prob = prob * probs[int(it-rvalues.begin())];
		}
	}
	

	return(prob);
}


// This  function is used for the proposal portion associated with a new
// sequence. For each locus it determines the PMF. Then either a new value
// for that locus is sampled based on the PMF and the probability is returned
// or the probability of the actual value, with respect to the PMF is returned.
double updateS(bool sample, treeNode* node_ptr,
				  const vector<int>& first_markers, double theta,
	  			  const vector<vector<double> >& allele_probs,
				  const vector<vector<double> >& hap_probs)
{
	int locus=0, locus_min=0, s_index=0, locus_max=0, j=0;
	double prob=1;
	char temp=' ';
	vector<double> probs;

	// Update the left side
	locus_min = first_markers[0]-node_ptr->z_L+1;
	for ( locus = first_markers[0]; locus > locus_min; locus-- ){

		getSijPMF("left", node_ptr, probs, first_markers, locus, theta,
					 allele_probs, hap_probs );


		// We are either sampling a new value for this locus based on Q() 
		// (sample=1) and must also return the prob. Or we are computing the
		// probability of the reverse proposal (ie finding Q() for the allele
		// that we have at this locus
		if ( sample == 1 ){ 
			s_index = samplePMF( probs ); 
			changeSij( node_ptr, itoc( s_index ), locus );
			prob = prob * probs[s_index];
			if (prob==0){
				/*cout*/ Rcpp::Rcout << "Prob = 0; Left; sample" << endl;
			}
		} else {
			temp = node_ptr->seq[locus];
			if (temp == '0'){
				prob = prob * probs[0];
			} else if (temp == '1'){
				prob = prob * probs[1];
			} else { /*cout*/ Rcpp::Rcout << "Neither 0 nor 1" << endl; }
			if (prob==0){
				/*cout*/ Rcpp::Rcout << "Prob = 0; Left; Prob" << endl;
			}
		}

	}

	if ( sample == 1 ){
		for (j=locus; j>=0; j--){
			node_ptr->seq[j]='-';
		}
	}

	// Update the right side
	locus_max = first_markers[1]+node_ptr->z_R-1;
	for ( locus = first_markers[1]; locus < locus_max; locus++ ){

		getSijPMF("right", node_ptr, probs, first_markers, locus, theta,
				  allele_probs, hap_probs );

		// We are either sampling a new value for this locus based on Q()
		// (sample=1) and must return the prob. Or we are computing the
		// probability of the reverse proposal (ie finding Q() for the allele
		// that we have at this locus
		if ( sample == 1 ){
			s_index = samplePMF( probs );
			changeSij( node_ptr, itoc( s_index ), locus );
			prob = prob * probs[s_index];
			if (prob==0){
				/*cout*/ Rcpp::Rcout << "Prob = 0; Right; sample" << endl;
			}
		} else {
			temp = node_ptr->seq[locus];
			if (temp == '0'){
				prob = prob * probs[0];
			} else if (temp == '1'){
				prob = prob * probs[1];
			} else { /*cout*/ Rcpp::Rcout << "Neither 0 nor 1" << endl; }
			if (prob==0){
				/*cout*/ Rcpp::Rcout << "Prob = 0; Right; prob" << endl;
			}
		}
		
	}

	if ( sample == 1 ){
		for (j=locus; j<(int)node_ptr->seq.size(); j++){
			node_ptr->seq[j]='-';
		}
	}
	
	return(prob);
	
}
	

// This is the modified proposal distribution for R
void getRiPMF( string side, treeNode* node_ptr, vector<double>& probs,
				vector<int>& validr, double rho,
				const vector<int>& first_markers, const vector<double>& locations,
				int ref, int type, treeNode* other_ptr,
				treeNode* sib_ptr)
{

	int r_sib=0, z_sib=0, r=0, z=0, r_a=0, z_a=0;
	int r_old_sib=0, /*z_old_sib=0,*/ r_other=0, z_other=0, z_gp=0;
	int init=1, end=0, i=0;
	double prob=1;


	// Resize probs to account for the support of r_i
	probs.clear();
	validr.clear();

	if (node_ptr->parent==NULL){
		probs.push_back(1);
		validr.push_back(0);
	}
	else {
		// Set up for the left and right side. Different variables are
		// accessed depending on side.
		if (side == "left"){
	
			z = node_ptr->z_L;
			r = node_ptr->r_L;
			r_a = node_ptr->parent->r_L;
			z_a = node_ptr->parent->z_L;
			if (sib_ptr != NULL){
				r_old_sib = sib_ptr->r_L;
//				z_old_sib = sib_ptr->z_L;
			}
			if (other_ptr != NULL){
				r_other = other_ptr->r_L;
				z_other = other_ptr->z_L;
			}
			if  ( (type == 5) || (type == 7) ){
				z_gp = node_ptr->parent->parent->z_L;
			}
			if (node_ptr->parent->child1->id == node_ptr->id){
				r_sib = node_ptr->parent->child2->r_L;
				z_sib = node_ptr->parent->child2->z_L;
			} else {
				r_sib = node_ptr->parent->child1->r_L;
				z_sib = node_ptr->parent->child1->z_L;
			}
		} else {
		
			z = node_ptr->z_R;
			r = node_ptr->r_R;
			r_a = node_ptr->parent->r_R;
			z_a = node_ptr->parent->z_R;
			if (sib_ptr != NULL){
				r_old_sib = sib_ptr->r_R;
//				z_old_sib = sib_ptr->z_R;
			}
			if (other_ptr != NULL){
				r_other = other_ptr->r_R;
				z_other = other_ptr->z_R;
			}
			if ( (type == 5) || (type == 7) ){
				z_gp = node_ptr->parent->parent->z_R;
			}
			if (node_ptr->parent->child1->id == node_ptr->id){
				r_sib = node_ptr->parent->child2->r_R;
				z_sib = node_ptr->parent->child2->z_R;
			} else {
				r_sib = node_ptr->parent->child1->r_R;
				z_sib = node_ptr->parent->child1->z_R;
			}
		}


		// Set up the possible values for r for each of the different types of
		// updates. These correspond to the local updates, the minor and major
		// topology rearrangements.
		end = z;
	
		switch (type) {
		case 1:
			if (z_sib<r_a){
				init = r_a;
			}
			break;
		case 2:
			if (r_sib < r_a){
				init = r_a;
			}
			break;
		case 3:
			if ( r_sib == z_a ){
				end = min(z_a,z);
			} else {
				init = z_a;
				end = z_a;
			}
			break;
		case 4:
			break;
		case 5:
			if ( (r_other != z_gp)&&(r_sib < z_gp ) ) {
				init = z_gp;	
			}
			break;
		case 6:
			if ( max(r_old_sib, z_other) >= z_a ){
				end = min(z_a,z);
			} else {
				init = z_a;
				end = z_a;
			}
			break;
		case 7:
			if ( (r_sib < z_gp) && (r_other < z_gp) ) {
				init = z_gp;
			}
			break;
		default:
			/*cerr*/ Rcpp::Rcout << "getRiPMF: Not a valid case" << endl;
			throw Rcpp::exception("getRiPMF: Not a valid case"); //exit(1);
		}


		// For each value of r_i in the support, compute the
		// numerator of the distribution
	
		for (i=init; i<=end; i++){
			if (side == "left"){
				node_ptr->r_L = i;
				prob = probRsided("left",node_ptr,rho,first_markers,locations,ref);
				probs.push_back(prob);
				validr.push_back(i);
			} else {
				node_ptr->r_R = i;
				prob = probRsided("right",node_ptr,rho,first_markers,locations,ref);
				probs.push_back(prob);
				validr.push_back(i);
			}
		}

		if ( side == "left" ){
			node_ptr->r_L = r;
		} else { node_ptr->r_R = r; }

	
	double sum=0;
	for (i=0; i<(int)probs.size(); i++){
		sum = sum + probs[i];
	}
	for (i=0; i<(int)probs.size(); i++){
		probs[i] = probs[i]/sum;
	}
	
	}
}



// This finds the probability distribution for an s_j^i
void getSijPMF(string side, treeNode* node_ptr, vector<double>& probs,
			   const vector<int>& first_markers, int locus, double theta,
	  		   const vector<vector<double> >& allele_probs,
   			   const vector<vector<double> >& hap_probs)
{

	int sij=0;
	unsigned int i=0;
	string seq = node_ptr->seq;
	double probA=1, probC1=1, probC2=1, sum=0;
	int rC1=0, rC2=0;
	
	probs.clear();

	// Find the index (with respect to the sequence) for the recombination
	// variables for both of the children of node_ptr
	if ( side == "left" ){
		rC1 = first_markers[0]-node_ptr->child1->r_L+1;
		rC2 = first_markers[0]-node_ptr->child2->r_L+1;
	} else {
		rC1 = first_markers[1]+node_ptr->child1->r_R-1;
		rC2 = first_markers[1]+node_ptr->child2->r_R-1;
	}
	
	
	for (sij=0; sij<2; sij++){

		// Change the locus to either 0 or 1
		node_ptr->seq[locus]=itoc(sij);

		// Find probability associated with similarity to the parent
		probA = probSij(side, node_ptr, locus, theta, first_markers, allele_probs,
						hap_probs);

		// Find probability associated with similarity to c1
		if ( ((side=="left")&&(locus>rC1)) || ((side=="right")&&(locus<rC1)) ){
			probC1 = probSij(side, node_ptr->child1, locus, theta, first_markers,
							 allele_probs, hap_probs);
		}
		
		// Find probability associated with similarity to c2
		if ( ((side=="left")&&(locus>rC2)) || ((side=="right")&&(locus<rC2)) ) {
			probC2=probSij(side, node_ptr->child2, locus, theta, first_markers,
						   allele_probs, hap_probs);
		}

		probs.push_back(probA*probC1*probC2);
	}
	
	node_ptr->seq = seq;

	for (i=0; i<probs.size(); i++){
		sum = sum + probs[i];
	}
	for (i=0; i<probs.size(); i++){
		probs[i] = probs[i]/sum;
	}

}


void updateI_rho( unsigned int& I, vector<double> constants,
		vector<double> lambdas, double& rho, //stringstream& ssOutFile, //ofstream& outfile,
		vector<treeNode*> node_ptrs, vector<int> first_markers,
		vector<double> locations, int ref_marker_location){

	unsigned int J=0;
	double qforward=1, qbackward=1, MHratio=0;
	double priorRatio=0, qRatio=0, PrRatio=0;
	double myunif = generateStdUnif();
//	int accept=0;

	if (lambdas[I-1]!=rho){
		/*cerr*/ Rcpp::Rcout << "There is a problem in updateI" << endl;
		throw Rcpp::exception("There is a problem in updateI"); //exit(1);
	}

	// Propose a new value for I. If I is at maximum, subtract 1. If I is at
	// minimum, add 1. Otherwise randomly choose to add or subtract 1. Note
	// that generateDiscreteUnif samples either x=0 or 1 here. If x=0, 2x-1=-1;
	// if x=1, 2x-1=1.
	if ( I==constants.size() ){
		J=I-1;
		qforward=1;
	} else if ( I==1 ){
		J=I+1;
		qforward=1;
	} else {
		J=I+(2*generateDiscreteUnif(0,1)-1);
		qforward=0.5;
	}

	if ( J==constants.size() ){
		qbackward=1;
	} else if ( J==1 ){
		qbackward=1;
	} else {
		qbackward=0.5;
	}

	qRatio = qbackward/qforward;

//	ssOutFile << I << " " << J << " "; //outfile

	for (int i=0; i<(int)node_ptrs.size()-1; i++){
		PrRatio = PrRatio +
				log(probRi( node_ptrs[i], lambdas[J-1], first_markers, locations,
					ref_marker_location )) -
				log(probRi( node_ptrs[i], lambdas[I-1], first_markers, locations,
					ref_marker_location ));
	}

	priorRatio = constants[J-1]/constants[I-1];


	 MHratio = min( exp(PrRatio) * priorRatio* qRatio , double( 1 ) );

	 if ( myunif < MHratio ){
		I=J;
		rho=lambdas[J-1];
//	 	accept = 1;
	 }

	 //outfile << qRatio << " " << PrAgivenGJ << " " << PrAgivenGI << " "
	 //		 << priorRatio << " ";
//	 ssOutFile << accept << " " << rho << endl; //outfile 

}


// This function computes the log of the posterior probability of all variables.

double probAll_relative(vector<treeNode*> node_ptrs, const double rho,
		const double theta, const vector<int>& first_markers,
		const vector<double>& locations, int ref_marker_location,
		const vector<vector<double> >& allele_probs,
		const vector<vector<double> >& hap_probs, // const int STflag,
		const double shp, const double scl )
{
	double logprob=0;
	int total_seq=(node_ptrs.size()+1)/2;
	int max_index=node_ptrs.size()-1;

	// Mutation rate model is uniform so don't return it. It just
	// adds a constant to all logprob.

	// If no ST, recombination rate has gamma distribution.

		logprob=logprob+ log(probGamma(rho, shp, scl));

	// Times are exponential
	logprob=logprob+ probT(node_ptrs, total_seq, max_index);

	// Traverse each node in nodeptrs and compute the probabilities of the R and S values
	for (int i=0; i<=max_index; i++){

		if (node_ptrs[i]->parent!=NULL){
			logprob=logprob+ log(probRi( node_ptrs[i], rho, first_markers, locations, ref_marker_location));
		}
		logprob=logprob+ log(probSi( node_ptrs[i], theta, first_markers, allele_probs, hap_probs));
	}
	return(logprob);
}




