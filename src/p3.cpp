/*----------------------------------------------------------------------------

  p3.cpp - This file contains the functions to complete the local updates.

  Kelly Burkett; July 2009

  Modification Notes:
  * Added values[2] = exp(values[2]); and values[3] = exp(values[3]); since
    the error check on the MH ratio assumes that these values have not had
    log's taken (Nov 2009)

  --------------------------------------------------------------------------*/


#include <iostream> //for cerr/cout
#include <cmath>
#include <Rcpp.h>

#include "nodefunctions.h"
#include "proposal.h"
#include "models.h"
#include "rng.h"
#include "treeBuild.h"
#include "p3.h"
#include "misc.h"


// This function does the local updates at each internal node in turn 
int updateP3( vector<treeNode*>& node_ptrs, vector<int> first_markers,
			   int ref_marker_location,
			   vector<double> locations, double rho, double theta,
			   vector<vector<double> > allele_probs,
			   vector<vector<double> > hap_probs, int it )
{
	unsigned int i=0;
	int accept = 0;
	treeNode *childptr1 = NULL, *childptr2 = NULL;
	double probA=1, probAnew=1, qA=1, qAnew=1, MHratio=1;
 	double myunif = generateStdUnif();
	string old_s = "";


	// Make a copy of the node_ptrs which has the update order (node_ptrs
	// itself will change order after an update is made because the order
	// of coalescent events may change
	vector<treeNode*> node_ptrs_order(node_ptrs.size());
	node_ptrs_order = node_ptrs;


	// Update each internal node by going to each ptr in node_ptrs rather than
	// by recursively traversing the tree. This will require that this vector
	// be kept in the order of coalescence events but the reordering must occur
	// after all times had been changed (otherwise, we could change the same
	// node twice.
	for ( i = (node_ptrs_order.size()+1)/2; i < node_ptrs_order.size(); i++ )
	{

		accept++;

		// Decide which child is updated first; it will be the oldest
		if ( node_ptrs_order[i]->child1->t > node_ptrs_order[i]->child2->t ){
			childptr1 = node_ptrs_order[i]->child1;
			childptr2 = node_ptrs_order[i]->child2;
		} else {
			childptr1 = node_ptrs_order[i]->child2;
			childptr2 = node_ptrs_order[i]->child1;
		}

		// Make copies of the nodes that are changed by this update since the
		// changes may be rejected. Here I am copying the value of treeNode
		// that is pointed to in to a new variable of type treeNode.
		treeNode copy_child1 = *(node_ptrs_order[i]->child1);
		treeNode copy_child2 = *(node_ptrs_order[i]->child2);
		treeNode copy_nodei = *node_ptrs_order[i];

		
		// Compute the term proportional to Pr(A|G) for this update
		probA  = probP3( node_ptrs_order[i], theta, rho, first_markers,
				 ref_marker_location, locations, allele_probs,
				 hap_probs, node_ptrs );

		// Computing terms in Q( A | tilde(A) )
		qA = qA * updateS( 0, node_ptrs_order[i], first_markers, theta,
						   allele_probs, hap_probs);
		qA = qA * updateR( 0, childptr1, rho, first_markers,
								 locations, ref_marker_location, 1 );
		qA = qA * updateR( 0, childptr2, rho, first_markers,
								 locations, ref_marker_location, 2 );
		qA = qA * updateT( 0, node_ptrs_order[i] );

								
		// Propose new t_i, r_c1, r_c2 and s^i; find Q( tilde(A) | A )
		qAnew = qAnew * updateT( 1, node_ptrs_order[i] );
		qAnew = qAnew * updateR( 1, childptr1, rho, first_markers,
								  locations, ref_marker_location, 1 );
		qAnew = qAnew * updateR( 1, childptr2, rho, first_markers,
								  locations, ref_marker_location, 2 );
		qAnew = qAnew * updateS( 1, node_ptrs_order[i], first_markers, theta,
						   allele_probs, hap_probs);

		
		// The order of coalescent events may have changed so reorder and
		// compute the term proportional to Pr(A~|G). The lower bound on the
		// change is the older of the two children, so min_index is 1 higher
		// than the index of that child. The upper bound is the index of the
		// parent.
		int min_index=max(childptr1->id,childptr2->id);
		int max_index=node_ptrs_order.size();//to deal with case of head node
						     //no sorting needed

		if (node_ptrs_order[i]->parent!=NULL){
			max_index = node_ptrs_order[i]->parent->id-1;
		}
		sortNodePtrs(node_ptrs, min_index, max_index);
		probAnew = probP3( node_ptrs_order[i], theta, rho, first_markers,
						   ref_marker_location, locations, allele_probs,
						   hap_probs, node_ptrs );

		
		// Check that there are no infinite or 0 values as this will cause an
		// undefined MH ratio
		bool errorFlag=0;

		errorFlag = checkValue( qAnew, 1 );
		if (errorFlag==1) {
			/*cout*/ Rcpp::Rcout << "WARNING: Probability of proposed local updates is infinity" << endl;
		}

		errorFlag = checkValue( qA, 0 );
		if (errorFlag==1) {
			Rcpp::Rcout << "ERROR: Probability of proposing current values in local updates is 0" << endl; // cerr << "ERROR: Probability of proposing current values in local updates is 0" << endl;
			throw Rcpp::exception("ERROR: Probability of proposing current values in local updates is 0"); //exit(1);
		}

		errorFlag = checkValue( exp(probAnew-probA), 1 );
		if (errorFlag==1) {
				/*cout*/ Rcpp::Rcout << "WARNING: Probability of data given proposed changes in local updates"
					 << "is infinity. Watch for getting stuck" << endl;
		}



		// Compute the MH ratio
		MHratio = min( exp(probAnew-probA+log(qA)-log(qAnew)), double( 1 ) );

		// Convert back to the old time configuration since we are not
		// accepting the change
		if ( myunif > MHratio ){

			changeR(node_ptrs_order[i]->child1, copy_child1);
			changeR(node_ptrs_order[i]->child2, copy_child2);
			changeT(node_ptrs_order[i], copy_nodei.t);
			changeS(node_ptrs_order[i], copy_nodei.seq );
			sortNodePtrs(node_ptrs, min_index, max_index);
			accept--;
		}

		
		// Reset all the values for the next iteration 
		qA=1;
		qAnew=1;
		probA=1;
		probAnew=1;
	
	}
	
	return( accept );  
	
}



// This function computes the term proportional to the log of the likelihood 
// for computing the acceptance probability of the update
double probP3( treeNode* node_ptr, double theta, double rho,
			   vector<int>& first_markers, int ref_marker_location,
	  		   vector<double>& locations, vector<vector<double> >& allele_probs,
			   vector<vector<double> >& hap_probs,
			   vector<treeNode*> node_ptrs )
{
	double prob = log(double(1));
	int i=0, min_index=0, max_index=0;

	treeNode* nodes[] = { node_ptr, node_ptr->child1, node_ptr->child2 };
	 
	for (i=0; i<3; i++){
		prob = prob + log(probSi( nodes[i], theta, first_markers, allele_probs,
							 hap_probs));
		prob = prob + log(probRi( nodes[i], rho, first_markers, locations,
							 ref_marker_location ));
	}
	
	if (node_ptr->child1->t == 0){//Offspring are terminal nodes
		min_index = (node_ptrs.size()+1)/2;
	} else { min_index = node_ptr->child1->id; }

	if (node_ptr->parent == NULL){//At head node
		max_index = node_ptr->id-1;
	} else { max_index = node_ptr->parent->id-1; }

	// NOTE: probT returns the log of the required value
	prob = prob + probT(node_ptrs, min_index, max_index );

	return( prob );

}


int updateP3_allnodes( vector<treeNode*>& node_ptrs, vector<int> first_markers,
			   int ref_marker_location,
			   vector<double> locations, double rho, double theta,
			   vector<vector<double> > allele_probs,
			   vector<vector<double> > hap_probs, int it )
{
	int i=0;
	int accept = 0;
	double probA=1, probAnew=1, qA=1, qAnew=1, MHratio=1;
 	double myunif = generateStdUnif();
	string old_s = "";
	int last=node_ptrs.size(), first=0;


	// Make a copy of the tree in case the update is not accepted
	treeNode* copy_head_node = copyTree(node_ptrs[node_ptrs.size()-1]);
	vector<treeNode*> copy_node_ptrs(node_ptrs.size());
	getNodePtrs(copy_node_ptrs, copy_head_node);


	// Update the times at all nodes.
	int num_seq = (node_ptrs.size()+1)/2;
	double time=0, totalT=0, rate=0;
	//double temp=0;

	for (i=0; i<num_seq-1; i++){

		// While there are num_seq-i sequences left determine the exponential
		// parameter and generate the times from it.
		rate = (num_seq-i) * (num_seq-i-1) / 2;
		time = generateExp(rate);
		//temp=rate*exp(-rate*time);
		//qAnew = qAnew*rate*exp(-rate*time);
		//cout << i << " time new " << temp << endl;

		totalT=totalT+time;
		changeT( node_ptrs[num_seq+i], totalT );
		//temp=rate*exp(-rate*(copy_node_ptrs[num_seq+i]->t-copy_node_ptrs[num_seq+i-1]->t));
		//qA = qA*rate*exp(-rate*(copy_node_ptrs[num_seq+i]->t-copy_node_ptrs[num_seq+i-1]->t));
		//cout << i << " time old " << temp << endl;

	}


	// Update each internal node by going to each ptr in node_ptrs rather than
	// by recursively traversing the tree. This will require that this vector
	// be kept in the order of coalescence events but the reordering must occur
	// after all times had been changed (otherwise, we could change the same
	// node twice.
	first=num_seq;

	for ( i = first; i < last; i++ )
	{

			qA = qA * updateR( 0, copy_node_ptrs[i]->child1, rho, first_markers,
					locations, ref_marker_location, 4 );
			qA = qA * updateR( 0, copy_node_ptrs[i]->child2, rho, first_markers,
					locations, ref_marker_location, 4 );
			qAnew = qAnew * updateR( 1, node_ptrs[i]->child1, rho, first_markers,
					locations, ref_marker_location, 4 );
			qAnew = qAnew * updateR( 1, node_ptrs[i]->child2, rho, first_markers,
					locations, ref_marker_location, 4 );
			qA = qA * updateS_child( 0, copy_node_ptrs[i], first_markers, theta,
					allele_probs, hap_probs);
			qAnew = qAnew * updateS_child( 1, node_ptrs[i], first_markers, theta,
					allele_probs, hap_probs);


	}

	probA= probP3_all( copy_node_ptrs, theta, rho, first_markers,
			ref_marker_location, locations, allele_probs,
			hap_probs );

	probAnew = probP3_all( node_ptrs, theta, rho, first_markers,
			ref_marker_location, locations, allele_probs,
			hap_probs);


	// Check that there are no infinite or 0 values as this will cause an
	// undefined MH ratio
	bool errorFlag=0;

	errorFlag = checkValue( qAnew, 1 );
	if (errorFlag==1) {
		/*cout*/ Rcpp::Rcout << "WARNING: Probability of proposed local updates is infinity" << endl;
	}

	errorFlag = checkValue( qA, 0 );
	if (errorFlag==1) {
		Rcpp::Rcout << "ERROR: Probability of proposing current values in local updates is 0" << endl; // cerr << "ERROR: Probability of proposing current values in local updates is 0" << endl;
		throw Rcpp::exception("ERROR: Probability of proposing current values in local updates is 0"); //exit(1);
	}


	errorFlag = checkValue( exp(probAnew-probA), 1 );
	if (errorFlag==1) {
		/*cout*/ Rcpp::Rcout << "WARNING: Probability of data given proposed changes in local updates"
				<< "is infinity. Watch for getting stuck" << endl;
	}


	// Compute the MH ratio
	MHratio = min( exp(probAnew-probA)*qA/qAnew, double( 1 ) );

	//cout << qA << " " << qAnew << " " << probA << " " << probAnew << " " << MHratio << endl;

	// Convert back to the old time configuration since we are not
	// accepting the change
	if ( myunif > MHratio ){

		deleteTree( node_ptrs[node_ptrs.size()-1] ) ;
		node_ptrs.clear();
		node_ptrs = copy_node_ptrs;
		accept = 0;

	} else {

		deleteTree( copy_node_ptrs[copy_node_ptrs.size()-1] );
		copy_node_ptrs.clear();
		accept=1;

	}



	return( accept );

}



// This function computes the term proportional to the log of the likelihood
// for computing the acceptance probability of the update
double probP3_all( vector<treeNode*> node_ptrs, double theta, double rho,
		vector<int>& first_markers, int ref_marker_location,
		vector<double>& locations, vector<vector<double> >& allele_probs,
		vector<vector<double> >& hap_probs )
{
	double prob = log(double(1));
	int i=0;
	int num_nodes=node_ptrs.size();


	for (i=0; i<num_nodes; i++){
		prob = prob + log(probSi( node_ptrs[i], theta, first_markers, allele_probs,
							 hap_probs));
		prob = prob + log(probRi( node_ptrs[i], rho, first_markers, locations,
							 ref_marker_location ));
	}


	// NOTE: probT returns the log of the required value
	//prob = prob + probT(node_ptrs, (num_nodes+1)/2, num_nodes-1 );

	return( prob );

}


// This  function is used for the proposal portion associated with a new
// sequence. For each locus it determines the PMF. Then either a new value
// for that locus is sampled based on the PMF and the probability is returned
// or the probability of the actual value, with respect to the PMF is returned.
double updateS_child(bool sample, treeNode* node_ptr,
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

		getSijPMF_child("left", node_ptr, probs, first_markers, locus, theta,
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

		getSijPMF_child("right", node_ptr, probs, first_markers, locus, theta,
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


// This finds the probability distribution for an s_j^i
void getSijPMF_child(string side, treeNode* node_ptr, vector<double>& probs,
			   const vector<int>& first_markers, int locus, double theta,
	  		   const vector<vector<double> >& allele_probs,
   			   const vector<vector<double> >& hap_probs)
{

	int sij=0;
	unsigned int i=0;
	string seq = node_ptr->seq;
	double probC1=1, probC2=1, sum=0;
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

		if ((probC1==1)&&(probC2==1)){
			Rcpp::Rcout << "I think there is a problem" << endl; // cerr << "I think there is a problem" << endl;
			throw Rcpp::exception("I think there is a problem"); //exit(1);
		}

		probs.push_back(probC1*probC2);
	}

	node_ptr->seq = seq;

	for (i=0; i<probs.size(); i++){
		sum = sum + probs[i];
	}
	for (i=0; i<probs.size(); i++){
		probs[i] = probs[i]/sum;
	}

}
