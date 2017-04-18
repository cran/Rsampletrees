/*----------------------------------------------------------------------------

  p4.cpp - This file contains functions needed to complete the major
	       topology rearrangement

  Kelly Burkett; July 2009

  Modification Notes:
  * Added values[2] = exp(values[2]); and values[3] = exp(values[3]); since
    the error check on the MH ratio assumes that these values have not had
    log's taken (Nov 2009)

  --------------------------------------------------------------------------*/


#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <Rcpp.h>

#include "treeBuild.h"
#include "nodefunctions.h"
#include "proposal.h"
#include "models.h"
#include "rng.h"
#include "p4.h"


// This function completes the major topology rearrangement (P4) by first
// determing the set of nodes that can be moved and sampling from that set.
// A second node that is compatible with the first is chosen and the 
// topology rearrangement is made. Other variables are changed as needed 
// and the acceptance probability is computed. If the proposal is rejected
// the old arrangement is restored by replacing the tree with a copy of the
// old tree. 
int updateP4( vector<treeNode*>& node_ptrs, double theta, double rho,
			  vector<int>& first_markers, int ref_marker,
	 		  vector<double>& locations, vector<vector<double> >& hap_probs,
	  		  vector<vector<double> >& allele_probs )
{

	int index=-1, accept = 1;
	double qA = 1, qAnew=1, MHratio = 1, probA = 1, probAnew = 1;
	double myunif = generateStdUnif();


	// Store the old values in case we switch back. Four nodes in total
	// are changed but it is easier to make a copy of the
	// entire tree. NOTE: I may want to consider a more efficient way in
	// the future
	treeNode* copy_head_node = copyTree(node_ptrs[node_ptrs.size()-1]);
	vector<treeNode*> copy_node_ptrs(node_ptrs.size());
	getNodePtrs(copy_node_ptrs, copy_head_node);


	// Find the set of nodes which can be sampled from and choose a node from
	// that set
	vector<treeNode*> valid_nodes(0);

	findValidNodesP4(valid_nodes, node_ptrs[node_ptrs.size()-1]);
	if ( valid_nodes.size()==0 ){
		// This update can't be done so return
		return( -1 );
	}
	index = generateDiscreteUnif(0,valid_nodes.size()-1);
	treeNode* node_ptr = valid_nodes[index];
	qAnew = qAnew * 1/double(valid_nodes.size());


	// Find the relatives of the chosen node
	treeNode* sib_ptr = NULL;
	treeNode* aunt_ptr = NULL;

	if (node_ptr->parent->child1->id == node_ptr->id){
		sib_ptr = node_ptr->parent->child2;
	} else {
		sib_ptr = node_ptr->parent->child1;
	}

	if (node_ptr->parent->parent->child1->id == node_ptr->parent->id){
		aunt_ptr = node_ptr->parent->parent->child2;
	} else {
		aunt_ptr = node_ptr->parent->parent->child1;
	}


	// Given the first node chosen, randomly select a second node based on its
	// similarity to the first one
	vector<double> probs(0);

	selectOtherNodeP4(probs, valid_nodes, node_ptrs, node_ptr, first_markers);

	if ( probs.size()==0 ){
		// This update can't be done as there is no other node to coalesce with
		// so return
		return( -1 );
	}
	index = samplePMF( probs );
	treeNode* other_ptr = valid_nodes[index];
	qAnew = qAnew * probs[index];


	// Find the sibling of the other node chosen
	treeNode* other_sib_ptr = NULL;

  	if (other_ptr->parent->child1->id == other_ptr->id){
		other_sib_ptr = other_ptr->parent->child2;
	} else {
		other_sib_ptr = other_ptr->parent->child1;
	}


	// Compute Pr(A|G)
	probA = probP4( node_ptr, other_ptr, sib_ptr, node_ptrs, theta, rho,
					first_markers, ref_marker, locations, allele_probs,
	 			    hap_probs );


	// Compute the terms in Q( A | tilde(A) ) corresponding to the t, r and s
	// variables
	qA = qA * updateT( 0, node_ptr->parent );

	qA = qA * updateR( 0, other_ptr,rho, first_markers, locations, ref_marker,
			   3 );
	qA = qA * updateR( 0, node_ptr, rho, first_markers, locations, ref_marker,
			   4 );
	qA = qA * updateR( 0, sib_ptr, rho, first_markers, locations, ref_marker,
			   5, aunt_ptr);
	qA = qA * updateR( 0, node_ptr->parent, rho, first_markers, locations,
			   ref_marker, 3);
	qA = qA * updateS( 0, node_ptr->parent, first_markers, theta,
			   allele_probs, hap_probs);


	// Change the topology of the tree as necessary and change the time
	// of the new node
	topologyChangeP4( node_ptr, other_ptr, sib_ptr, node_ptrs );


	// Generate new r's for K_j (sib_ptr), K_i (node_ptr), K_c (other_ptr),
	// and K_p (node_ptr->parent). Generate a new s for K_p
	qAnew = qAnew * updateT( 1, node_ptr->parent );

	qAnew = qAnew * updateR( 1, sib_ptr, rho, first_markers, locations,
							 ref_marker, 3 );
	qAnew = qAnew * updateR( 1, node_ptr, rho, first_markers, locations,
							 ref_marker, 4 );
	qAnew = qAnew * updateR( 1, other_ptr, rho, first_markers, locations,
							 ref_marker, 5, other_sib_ptr);
	qAnew = qAnew * updateR( 1, node_ptr->parent, rho, first_markers,
							 locations, ref_marker, 3);
	qAnew = qAnew * updateS( 1, node_ptr->parent, first_markers, theta,
					  		 allele_probs, hap_probs);



	// Since the times are changed, the order of coalescent events may have
	// changed. Must re-order the vector of node pointers
	int min_index = min( max (node_ptr->id,other_ptr->id),
						 max(node_ptr->id,sib_ptr->id) );
	int max_index = max( node_ptr->parent->parent->id, sib_ptr->parent->id )-1;

	sortNodePtrs(node_ptrs, min_index, max_index);


	// Find the probability for choosing node K_i based on tilde(A). This is
	// used to compute qA
	valid_nodes.clear();
	findValidNodesP4(valid_nodes, node_ptrs[node_ptrs.size()-1]);
	qA = qA * 1/double(valid_nodes.size());


	// Find the probability for choosing the other node K_j based on
	// tilde(A) which is used in computing Q( A|tilde(A) )

	selectOtherNodeP4(probs, valid_nodes, node_ptrs, node_ptr, first_markers);

	vector<treeNode*>::iterator it;
	it = find( valid_nodes.begin(), valid_nodes.end(), sib_ptr);
	if ( it == valid_nodes.end() ){
		Rcpp::Rcout << "ERROR: Reverse change is not possible in Major Topology Change" << endl; // cerr << "ERROR: Reverse change is not possible in Major" << "Topology Change" << endl;
		throw Rcpp::exception("ERROR: Reverse change is not possible in Major Topology Change"); //exit(1);
	}
	qA = qA * probs.at(int( it - valid_nodes.begin() ));




	// Compute Pr(tilde(A)|G)
	probAnew = probP4( node_ptr, sib_ptr, other_ptr, node_ptrs, theta, rho,
					   first_markers, ref_marker, locations, allele_probs,
					   hap_probs );


	// Check that there are no infinite or 0 values as this will cause an
	// undefined MH ratio
	bool errorFlag=0;

	errorFlag = checkValue( qAnew, 1 );
	if (errorFlag==1) {
		/*cout*/ Rcpp::Rcout << "WARNING: Probability of proposed change in major topology update is infinity" << endl;
	}

	errorFlag = checkValue( qA, 0 );
	if (errorFlag==1) {
		Rcpp::Rcout << "ERROR: Probability of proposing current values in major topology update is 0" << endl; // cerr << "ERROR: Probability of proposing current values in major topology update is 0" << endl;
		throw Rcpp::exception("ERROR: Probability of proposing current values in major topology update is 0"); //exit(1);
	}

	errorFlag = checkValue( exp(probAnew-probA), 1 );
	if (errorFlag==1) {
			/*cout*/ Rcpp::Rcout << "WARNING: Probability of data given proposed changes in major topology update"
				 << "is infinity. Watch for getting stuck" << endl;
	}


	// Compute the MH ratio for this update
	MHratio = min( exp(probAnew-probA+log(qA)-log(qAnew)), double( 1 ) );


	// Convert back to the old configuration since we are not
	// accepting the change. The topology change is fairly big, so it is
	// easiest to just have the set of node pointers point to the copy
	// that was made at the beginning of this update.
	if ( myunif > MHratio ){

		deleteTree( node_ptrs[node_ptrs.size()-1] ) ;
		node_ptrs.clear();
		node_ptrs = copy_node_ptrs;
		accept = 0;

	} else {

		deleteTree( copy_node_ptrs[copy_node_ptrs.size()-1] );
		copy_node_ptrs.clear();

	}


	return( accept );
}


int updateP4_2( vector<treeNode*>& node_ptrs, double theta, double rho,
			  vector<int>& first_markers, int ref_marker,
	 		  vector<double>& locations, vector<vector<double> >& hap_probs,
	  		  vector<vector<double> >& allele_probs, int reps)
{

	int index=-1, accept = 1;
	double qA = 1, qAnew=1, MHratio = 1, probA = 1, probAnew = 1;
	double myunif = generateStdUnif(), temp=0;
	bool outflag=0;

if ( (reps>500000)&&(reps<510000) ){
	outflag=1;
}
	
	// Store the old values in case we switch back. Four nodes in total
	// are changed but it is easier to make a copy of the
	// entire tree. NOTE: I may want to consider a more efficient way in
	// the future
	treeNode* copy_head_node = copyTree(node_ptrs[node_ptrs.size()-1]);
	vector<treeNode*> copy_node_ptrs(node_ptrs.size());
	getNodePtrs(copy_node_ptrs, copy_head_node);

	
	// Find the set of nodes which can be sampled from and choose a node from
	// that set
	vector<treeNode*> valid_nodes(0);

	findValidNodesP4(valid_nodes, node_ptrs[node_ptrs.size()-1]);
	if ( valid_nodes.size()==0 ){
		// This update can't be done so return
		return( -1 );
	}
	index = generateDiscreteUnif(0,valid_nodes.size()-1);
	treeNode* node_ptr = valid_nodes[index];
	qAnew = qAnew * 1/double(valid_nodes.size());

	
	// Find the relatives of the chosen node
	treeNode* sib_ptr = NULL;
	treeNode* aunt_ptr = NULL;
	
	if (node_ptr->parent->child1->id == node_ptr->id){
		sib_ptr = node_ptr->parent->child2;
	} else {
		sib_ptr = node_ptr->parent->child1;
	}

	if (node_ptr->parent->parent->child1->id == node_ptr->parent->id){
		aunt_ptr = node_ptr->parent->parent->child2;
	} else {
		aunt_ptr = node_ptr->parent->parent->child1;
	}

	
	// Given the first node chosen, randomly select a second node based on its
	// similarity to the first one
	vector<double> probs(0);
	
	selectOtherNodeP4_2(probs, valid_nodes, node_ptrs, node_ptr, first_markers,
			theta, rho, locations);
	if ( probs.size()==0 ){
		// This update can't be done as there is no other node to coalesce with
		// so return
		return( -1 );
	}
	index = samplePMF( probs );
	treeNode* other_ptr = valid_nodes[index];
	qAnew = qAnew * probs[index];

if (outflag==1){
	/*cout*/ Rcpp::Rcout << endl << endl;
	/*cout*/ Rcpp::Rcout << "*********************************************************" << endl;
	/*cout*/ Rcpp::Rcout << "Chosen nodes sequence: " << node_ptr->seq << endl;
	/*cout*/ Rcpp::Rcout << "Old sib sequence: " << sib_ptr->seq << endl;
	/*cout*/ Rcpp::Rcout << "Sequences and Probs" << endl;
	for (unsigned int k=0; k<probs.size(); k++){
		/*cout*/ Rcpp::Rcout << valid_nodes[k]->seq << " " << probs[k] << endl;
	}
	/*cout*/ Rcpp::Rcout << endl;
	/*cout*/ Rcpp::Rcout << "New sib sequence: " << other_ptr->seq << " " << probs[index] << endl;
}
	
	// Find the sibling of the other node chosen
	treeNode* other_sib_ptr = NULL;
	
  	if (other_ptr->parent->child1->id == other_ptr->id){
		other_sib_ptr = other_ptr->parent->child2;
	} else {
		other_sib_ptr = other_ptr->parent->child1;
	}

	
	// Compute Pr(A|G)
if (outflag==1){
 	/*cout*/ Rcpp::Rcout << "Pr(A|G): " << endl;
}
	probA = probP4( node_ptr, other_ptr, sib_ptr, node_ptrs, theta, rho,
					first_markers, ref_marker, locations, allele_probs,
	 			    hap_probs, outflag );


	
	// Compute the terms in Q( A | tilde(A) ) corresponding to the t, r and s
	// variables
if (outflag==1){
	/*cout*/ Rcpp::Rcout << "Q(A~->A): " << endl;
}

	qA = qA * updateT( 0, node_ptr->parent );
if (outflag==1){
	temp = updateT( 0, node_ptr->parent );
	/*cout*/ Rcpp::Rcout << "updateT: " << temp << endl;
}

	qA = qA * updateR( 0, other_ptr,rho, first_markers, locations, ref_marker,
			   3 );
if (outflag==1){
	temp = updateR( 0, other_ptr,rho, first_markers, locations, ref_marker,
				   3 );
	/*cout*/ Rcpp::Rcout << "updateR (new sib): " << temp << endl;
}

	qA = qA * updateR( 0, node_ptr, rho, first_markers, locations, ref_marker,
			   4 );
if (outflag==1){
	temp=updateR( 0, node_ptr, rho, first_markers, locations, ref_marker,
						   4 );
	/*cout*/ Rcpp::Rcout << "updateR (moved node): " << temp << endl;
}

	qA = qA * updateR( 0, sib_ptr, rho, first_markers, locations, ref_marker,
			   5, aunt_ptr);
if (outflag==1){
	temp=updateR( 0, sib_ptr, rho, first_markers, locations, ref_marker,
						   5, aunt_ptr);
	/*cout*/ Rcpp::Rcout << "updateR (old sib): " << temp << endl;
}

	qA = qA * updateR( 0, node_ptr->parent, rho, first_markers, locations,
			   ref_marker, 3);
if (outflag==1){
	temp=updateR( 0, node_ptr->parent, rho, first_markers, locations,
						   ref_marker, 3);
	/*cout*/ Rcpp::Rcout << "updateR (moved node parent): " << temp << endl;
}

	qA = qA * updateS( 0, node_ptr->parent, first_markers, theta,
			   allele_probs, hap_probs);
if (outflag==1){
	temp=updateS( 0, node_ptr->parent, first_markers, theta,
						   allele_probs, hap_probs);
	/*cout*/ Rcpp::Rcout << "updateS (moved node parent): " << temp << endl;
}


	// Change the topology of the tree as necessary and change the time
	// of the new node
	topologyChangeP4( node_ptr, other_ptr, sib_ptr, node_ptrs );


	// Generate new r's for K_j (sib_ptr), K_i (node_ptr), K_c (other_ptr),
	// and K_p (node_ptr->parent). Generate a new s for K_p
if (outflag==1){
		/*cout*/ Rcpp::Rcout << "Q(A->A~): " << endl;
}


	qAnew = qAnew * updateT( 1, node_ptr->parent );
if (outflag==1){
	temp=updateT( 1, node_ptr->parent );
	/*cout*/ Rcpp::Rcout << "updateT: " << temp << endl;
}

	qAnew = qAnew * updateR( 1, sib_ptr, rho, first_markers, locations,
							 ref_marker, 3 );
if (outflag==1){
	temp=updateR( 1, sib_ptr, rho, first_markers, locations,
			 ref_marker, 3 );
	/*cout*/ Rcpp::Rcout << "updateR (old sib): " << temp << endl;
}

	qAnew = qAnew * updateR( 1, node_ptr, rho, first_markers, locations,
							 ref_marker, 4 );
if (outflag==1){
	temp=updateR( 1, node_ptr, rho, first_markers, locations,
			 ref_marker, 4 );
	/*cout*/ Rcpp::Rcout << "updateR (moved node): " << temp << endl;
}

	qAnew = qAnew * updateR( 1, other_ptr, rho, first_markers, locations,
							 ref_marker, 5, other_sib_ptr);
if (outflag==1){
	temp=updateR( 1, other_ptr, rho, first_markers, locations,
			 ref_marker, 5, other_sib_ptr);
	/*cout*/ Rcpp::Rcout << "updateR (new sib): " << temp << endl;
}

	qAnew = qAnew * updateR( 1, node_ptr->parent, rho, first_markers,
							 locations, ref_marker, 3);
if (outflag==1){
	temp=updateR( 1, node_ptr->parent, rho, first_markers,
			 locations, ref_marker, 3);
	/*cout*/ Rcpp::Rcout << "updateR (moved node parent): " << temp << endl;
}

	qAnew = qAnew * updateS( 1, node_ptr->parent, first_markers, theta,
					  		 allele_probs, hap_probs);
if (outflag==1){
	temp=updateS( 1, node_ptr->parent, first_markers, theta,
	  		 allele_probs, hap_probs);
	/*cout*/ Rcpp::Rcout << "updateS (moved node parent): " << temp << endl;
}

	
	// Since the times are changed, the order of coalescent events may have
	// changed. Must re-order the vector of node pointers
	int min_index = min( max (node_ptr->id,other_ptr->id),
						 max(node_ptr->id,sib_ptr->id) );
	int max_index = max( node_ptr->parent->parent->id, sib_ptr->parent->id )-1;
	
	sortNodePtrs(node_ptrs, min_index, max_index);

	
	// Find the probability for choosing node K_i based on tilde(A). This is
	// used to compute qA
	valid_nodes.clear();
	findValidNodesP4(valid_nodes, node_ptrs[node_ptrs.size()-1]);
	qA = qA * 1/double(valid_nodes.size());

	
	// Find the probability for choosing the other node K_j based on
	// tilde(A) which is used in computing Q( A|tilde(A) )

	selectOtherNodeP4_2(probs, valid_nodes, node_ptrs, node_ptr, first_markers,
			theta, rho, locations);
	vector<treeNode*>::iterator it;
	it = find( valid_nodes.begin(), valid_nodes.end(), sib_ptr);
	if ( it == valid_nodes.end() ){
		Rcpp::Rcout << "ERROR: Reverse change is not possible in Major Topology Change" << endl; // cerr << "ERROR: Reverse change is not possible in Major Topology Change" << endl;
		throw Rcpp::exception("ERROR: Reverse change is not possible in Major Topology Change"); //exit(1);
	}
	qA = qA * probs.at(int( it - valid_nodes.begin() ));


if (outflag==1){
	/*cout*/ Rcpp::Rcout << "Q(A~->A) Old sib sequence prob: " <<  probs.at(int( it - valid_nodes.begin() )) << endl;
}

	// Compute Pr(tilde(A)|G)
if (outflag==1){
 	/*cout*/ Rcpp::Rcout << "Pr(A~|G): " << endl;
}
	probAnew = probP4( node_ptr, sib_ptr, other_ptr, node_ptrs, theta, rho,
					   first_markers, ref_marker, locations, allele_probs,
					   hap_probs, outflag );

	
	// Check that there are no infinite or 0 values as this will cause an
	// undefined MH ratio
	bool errorFlag=0;

	errorFlag = checkValue( qAnew, 1 );
	if (errorFlag==1) {
		/*cout*/ Rcpp::Rcout << "WARNING: Probability of proposed change in major topology update is infinity" << endl;
	}

	errorFlag = checkValue( qA, 0 );
	if (errorFlag==1) {
		Rcpp::Rcout << "ERROR: Probability of proposing current values in major topology update is 0" << endl; // cerr << "ERROR: Probability of proposing current values in major topology update is 0" << endl;
		throw Rcpp::exception("ERROR: Probability of proposing current values in major topology update is 0"); //exit(1);
	}

	errorFlag = checkValue( exp(probAnew-probA), 1 );
	if (errorFlag==1) {
			/*cout*/ Rcpp::Rcout << "WARNING: Probability of data given proposed changes in major topology update"
				 << "is infinity. Watch for getting stuck" << endl;
	}

	
	// Compute the MH ratio for this update
	MHratio = min( exp(probAnew-probA+log(qA)-log(qAnew)), double( 1 ) );

	
	// Convert back to the old configuration since we are not
	// accepting the change. The topology change is fairly big, so it is
	// easiest to just have the set of node pointers point to the copy
	// that was made at the beginning of this update.
	if ( myunif > MHratio ){

		deleteTree( node_ptrs[node_ptrs.size()-1] ) ;
		node_ptrs.clear();
		node_ptrs = copy_node_ptrs;
		accept = 0;

	} else {

		deleteTree( copy_node_ptrs[copy_node_ptrs.size()-1] );
		copy_node_ptrs.clear();

	}

if (outflag==1){
	/*cout*/ Rcpp::Rcout << "alpha: " << MHratio << "; accept: " << accept << endl;
	/*cout*/ Rcpp::Rcout << endl << endl;
}
	return( accept );
}
// This function is used to create a vector of node pointers that point to
// nodes that can moved from their current location without causing
// inconsistencies in the topology. These are nodes where the maximum length
// of sequence is carried to a node's grandparent through either its sibling
// or its aunt, for both the left and right hand sequence.
// NOTE: With node_ptrs, this could be replaced with a non-recursive function
// to make it easier to follow
void findValidNodesP4(vector<treeNode*>& valid_nodes, treeNode* head_node)
{
	treeNode* sib_ptr = NULL;
	treeNode* aunt_ptr = NULL;

	if (head_node != NULL){

		// The head node and both its children are not valid for moving
		if ((head_node->parent == NULL)||(head_node->parent->parent == NULL)){
			
			findValidNodesP4(valid_nodes, head_node->child1);
			findValidNodesP4(valid_nodes, head_node->child2);
			
		} else {

			// Find the relatives
			if (head_node->parent->child1->id == head_node->id){
				sib_ptr = head_node->parent->child2;
			} else {
				sib_ptr = head_node->parent->child1;
			}

			if (head_node->parent->parent->child1->id == head_node->parent->id){
				aunt_ptr = head_node->parent->parent->child2;
			} else {
				aunt_ptr = head_node->parent->parent->child1;
			}

			// Determine whether this node can be moved. 
			if (head_node->parent->parent->z_L <= max(aunt_ptr->r_L,sib_ptr->z_L)){
				if (head_node->parent->parent->z_R <= max(aunt_ptr->r_R,sib_ptr->z_R)){
					valid_nodes.push_back(head_node);
				}
			}
		
			findValidNodesP4(valid_nodes, head_node->child1);
			findValidNodesP4(valid_nodes, head_node->child2);
		}

	}

	
}


// This function determines a distribution for sampling another node
// based on sequence similarity to the chosen node, and stores those nodes
// in a vector called values. Only a node that is compatible for the topology
// rearrangement will be added to the set of nodes that can be sampled from
// (the time of the parent must be older than the chosen node, it can't be a
// parent, sib, or aunt of the chosen node)
void selectOtherNodeP4(vector<double>& probs, vector<treeNode*>& values,
					   vector<treeNode*>& node_ptrs, treeNode* chosen_node,
					   vector<int>& first_markers)
{

	unsigned int i=0;
	double total=0, sim=0, epsilon=0.0001;
	int l_min=0, r_min=0, l_index=0, r_index=0, j=0, good_node=1;
	
	probs.clear();
	values.clear();
	
	for (i=0; i<node_ptrs.size()-1; i++){

		// Determine if this node is eligible to be chosen. There are numerous
		// conditions that will exclude this node, hence the many if's. They
		// could all be combined but it was very confusing to figure out
		// what the condition was

		// Current node is the same as the chosen node
		if (node_ptrs[i]->id == chosen_node->id){
			good_node=0;
		}

		// Both nodes are sibs
		else if (node_ptrs[i]->parent->id == chosen_node->parent->id){
			good_node=0;
		}

		// Aunt/Niece relationship
		else if ( (node_ptrs[i]->parent->id == chosen_node->parent->parent->id) ||
			 ((node_ptrs[i]->parent->parent!=NULL)&&
			  (node_ptrs[i]->parent->parent->id == chosen_node->parent->id)) ){
			good_node=0;
		}

		// Current node is parent/grandparent of the chosen node (note that the
		// opposite relationship is not allowed by the condition on time)
		else if ( (node_ptrs[i]->id == chosen_node->parent->id) ||
			 (node_ptrs[i]->id == chosen_node->parent->parent->id) ){
			good_node=0;
		}

		// Time condition on the node. If the parent of the current node
		// is younger than the chosen node, it is not possible to put the
		// chosen node's parent less than current node's parent 
		else if (node_ptrs[i]->parent->t <= chosen_node->t) {
			good_node=0;
		}

		if (good_node == 1){
			
			sim = 0;
			values.push_back(node_ptrs[i]);
		
			/*if (node_ptrs[i]->parent==NULL){
				l_min = min(chosen_node->r_L, node_ptrs[i]->z_L)-1;
				r_min = min(chosen_node->r_R, node_ptrs[i]->z_R)-1;
			} else {*/
				l_min = min(chosen_node->z_L, node_ptrs[i]->z_L)-1;
				r_min = min(chosen_node->z_R, node_ptrs[i]->z_R)-1;
			//}


			if ( (l_min==0)&&(r_min==0) ){
				// If the sequence is not passed down to its ancestors at all due
				// to recombination events on either side, then there is no
				// sequence to compare, so let the probability of being chosen be
				// small (epsilon)
				probs.push_back(epsilon);

			} else {
				// Otherwise, we know the sequence at this node and can compare
				// the similarity. Note that if the two sequences are
				// completely different, the value for sim will also be epsilon
				l_index = first_markers[0]-l_min+1;
				r_index = first_markers[1]+r_min-1;

				for (j=l_index; j<=r_index; j++){
					if (node_ptrs[i]->seq[j] == chosen_node->seq[j]){
						sim+=1;
					}
				}

				if (sim==0){
					probs.push_back(epsilon);
				} else{
					probs.push_back( (sim+epsilon)/(l_min+r_min+epsilon) );
				}

			}
			total += probs[probs.size()-1];
				
		}

		good_node=1;

	}

	for (i=0; i<probs.size(); i++){
		probs[i] = probs[i]/total;
	}
	
			
	
}


void selectOtherNodeP4_2(vector<double>& probs, vector<treeNode*>& values,
					   const vector<treeNode*>& node_ptrs, treeNode* chosen_node,
					   const vector<int>& first_markers, double theta, double rho,
					   const vector<double>& locations)
{

	unsigned int i=0;
	double total=0, /*sim=0,*/ epsilon=0.000001, lambda=0, t=0, chosenl, chosenr, otherl, otherr, prob;
	int lmin=0, rmin=0, l_index=0, r_index=0, j=0, good_node=1;

	probs.clear();
	values.clear();

	for (i=0; i<node_ptrs.size()-1; i++){

		// Determine if this node is eligible to be chosen. There are numerous
		// conditions that will exclude this node, hence the many if's. They
		// could all be combined but it was very confusing to figure out
		// what the condition was

		// Current node is the same as the chosen node
		if (node_ptrs[i]->id == chosen_node->id){
			good_node=0;
		}

		// Both nodes are sibs
		else if (node_ptrs[i]->parent->id == chosen_node->parent->id){
			good_node=0;
		}

		// Aunt/Niece relationship
		else if ( (node_ptrs[i]->parent->id == chosen_node->parent->parent->id) ||
			 ((node_ptrs[i]->parent->parent!=NULL)&&
			  (node_ptrs[i]->parent->parent->id == chosen_node->parent->id)) ){
			good_node=0;
		}

		// Current node is parent/grandparent of the chosen node (note that the
		// opposite relationship is not allowed by the condition on time)
		else if ( (node_ptrs[i]->id == chosen_node->parent->id) ||
			 (node_ptrs[i]->id == chosen_node->parent->parent->id) ){
			good_node=0;
		}

		// Time condition on the node. If the parent of the current node
		// is younger than the chosen node, it is not possible to put the
		// chosen node's parent less than current node's parent
		else if (node_ptrs[i]->parent->t <= chosen_node->t) {
			good_node=0;
		}

		if (good_node == 1){

			values.push_back(node_ptrs[i]);

			// Estimate the branch length separating these two nodes
			t = (node_ptrs[i]->parent->t + max(chosen_node->t,node_ptrs[i]->t))/2;
			lambda=(2*t-chosen_node->t-node_ptrs[i]->t)*theta/2;

			// Estimate the location of recombination break points given the
			// estimated branch length. Use the expected value of the exp
			// distribution: 1/lambda

			chosenl=locations[first_markers[0]]-1/(rho/2*(t-chosen_node->t));
			if (chosenl<(locations[first_markers[0]-chosen_node->z_L+2])){
				chosenl=locations[first_markers[0]-chosen_node->z_L+2];
			}
			otherl=locations[first_markers[0]]-1/(rho/2*(t-node_ptrs[i]->t));
			if (otherl<(locations[first_markers[0]-node_ptrs[i]->z_L+2])){
				otherl=locations[first_markers[0]-node_ptrs[i]->z_L+2];
			}
			lmin= max(chosenl,otherl);
			l_index=first_markers[0];

			while ( (l_index>=0)&&(locations[l_index]>=lmin) ){
				l_index--;
			}
			l_index++;

			chosenr=locations[first_markers[1]]+1/(rho/2*(t-chosen_node->t));
			if (chosenr>(locations[first_markers[1]+chosen_node->z_R-2])){
				chosenr=locations[first_markers[1]+chosen_node->z_R-2];
			}
			otherr=locations[first_markers[1]]+1/(rho/2*(t-node_ptrs[i]->t));
			if (otherr>(locations[first_markers[1]+node_ptrs[i]->z_R-2])){
				otherr=locations[first_markers[1]+node_ptrs[i]->z_R-2];
			}
			rmin= min(chosenr,otherr);
			r_index=first_markers[1];


			while ( (r_index<(int)locations.size())&&(locations[r_index]<=rmin) ){
				r_index++;
			}
			r_index--;


			/*lmin = min(chosen_node->z_L, node_ptrs[i]->z_L)-1;
			rmin = min(chosen_node->z_R, node_ptrs[i]->z_R)-1;
			l_index = first_markers[0]-lmin+2;
			r_index = first_markers[1]+rmin-2;*/

			if ( (l_index==(first_markers[0]+1))&&(r_index==(first_markers[1]-1)) ){
			//if ( (lmin==1)&&(rmin==1) ){
				// If the sequence is not passed down to its ancestors at all due
				// to recombination events on either side, then there is no
				// sequence to compare, so let the probability of being chosen be
				// small (epsilon)
				probs.push_back(epsilon);

			} else {

				prob=1;
				for (j=l_index; j<=r_index; j++){
					if (node_ptrs[i]->seq[j] == chosen_node->seq[j]){
						prob = prob*(1/ALPHA*(1-exp(-lambda))+exp(-lambda));
					} else {
						prob = prob*1/ALPHA*(1-exp(-lambda));
					}
				}
				probs.push_back(prob);


			}
			total += probs[probs.size()-1];

		}

		good_node=1;

	}

	for (i=0; i<probs.size(); i++){
		probs[i] = probs[i]/total;
	}


}




// This function is used to rearrange the nodes for the major topology
// rearrangement.
void topologyChangeP4( treeNode* node_ptr, treeNode* chosen_ptr,
					   treeNode* sib_ptr, vector<treeNode*>& node_ptrs )
{
	
	// Change the ptrs of the sibling
	sib_ptr->parent = sib_ptr->parent->parent;
	if (node_ptr->parent->parent->child1->id == node_ptr->parent->id){
		sib_ptr->parent->child1 = sib_ptr;
	} else {
		sib_ptr->parent->child2 = sib_ptr;
	}
	
	// Change ptr of c's parents child to point to p rather than c
	if (chosen_ptr->parent->child1->id == chosen_ptr->id){
		chosen_ptr->parent->child1 = node_ptr->parent;
	} else {
		chosen_ptr->parent->child2 = node_ptr->parent;
	}

	// Make p's parent c's current parent
	node_ptr->parent->parent = chosen_ptr->parent;

	// Change c's parent to be p
	chosen_ptr->parent = node_ptr->parent;

	// Make c one of p's children
	if (node_ptr->parent->child1->id == node_ptr->id){
		node_ptr->parent->child2 = chosen_ptr;
	} else {
		node_ptr->parent->child1 = chosen_ptr;
	}
	
	// Update the branch length for sib_ptr because it has a new parent
	sib_ptr->b = sib_ptr->parent->t-sib_ptr->t;

}


// This function computes a term proportional to the  log of the likelikood
// of the augmented data. 
double probP4(treeNode* node_ptr, treeNode* other_ptr, treeNode* sib_ptr,
			  vector<treeNode*> node_ptrs,
	 		  double theta, double rho, vector<int>& first_markers,
	  		  int ref_marker, vector<double>& locations,
	  		  vector<vector<double> >& allele_probs, vector<vector<double> >& hap_probs,
	  		  bool outflag)
{
	double prob = log(double(1)), temp=0;

	// Compute the prob(omega) terms. For the major topology rearrangement,
	// the minimum intercoalescence time that could change is between the youngest
	// of the two sets of children and the next coalescent event. If this child is
	// K_min, then the intercoalescence time is t_{min+1}-t_min. Thus we pass min+1
	// to the function since the first term it computes will be t_{min+1}-t_min.
	// This is achieved by not subtracting 1 since min_index will index the vector
	// of node_ptrs.
	int min_index = min( max(node_ptr->id,other_ptr->id), max(node_ptr->id,sib_ptr->id) );

	// Since the ID's and the vector of node pointers are in order of coalescence,
	// the ID of the node-1 is used to index the vector of node pointers.
	int max_index = max( node_ptr->parent->parent->id, other_ptr->parent->id )-1;

if (outflag==1){
	temp =  probT( node_ptrs, min_index, max_index );
	/*cout*/ Rcpp::Rcout << "probT: " << exp(temp) << endl;
}
	prob = prob + probT( node_ptrs, min_index, max_index );
	
	// Compute the prob(S) and Prob(R) terms for nodes i (node_ptr), j (sib_ptr),
	// c (other_ptr) and p (node_ptr->parent)
	treeNode* nodes[] = {node_ptr, other_ptr, sib_ptr, node_ptr->parent};

	for (int i=0; i<4; i++){
		prob = prob + log(probSi( nodes[i], theta, first_markers, allele_probs,
				hap_probs ));
		prob = prob + log(probRi( nodes[i], rho, first_markers, locations,
							 ref_marker ));
if (outflag==1){
	temp =  probSi( nodes[i], theta, first_markers, allele_probs,
			hap_probs );
	/*cout*/ Rcpp::Rcout << "probS: " << temp << endl;
	temp = probRi( nodes[i], rho, first_markers, locations,
			 ref_marker );
}
	}

	return( prob );
	
}
