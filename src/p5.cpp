/*----------------------------------------------------------------------------

  p5.cpp - This file contains functions needed to complete the minor
		 topology rearrangement

  Kelly Burkett; July 2009

  Modification Notes:

  --------------------------------------------------------------------------*/

#include <vector>
#include <cmath>
#include <iostream>
#include <Rcpp.h>

#include "treeBuild.h"
#include "nodefunctions.h"
#include "proposal.h"
#include "models.h"
#include "rng.h"
#include "p5.h"



// This is the update function for completing the minor topology rearrangement
// (P5).
int updateP5( vector<treeNode*>& node_ptrs, double theta, double rho,
			  vector<int> first_markers, int ref_marker, vector<double> locations,
	 		  vector<vector<double> > hap_probs,
	  		  vector<vector<double> > allele_probs)
{
	double qA = 1, qAnew =1, MHratio = 1, probA = 1, probAnew = 1;
	int accept = 1;
	double myunif = generateStdUnif();

	
	// Find the set of nodes which can be sampled from and choose a node from
	// them.
	vector<treeNode*> valid_nodes(0);

	findValidNodesP5(valid_nodes, node_ptrs[node_ptrs.size()-1]);
	if ( valid_nodes.size()==0 ){
		// This update can't be done so return
		return( 0 );
	}
	int index = generateDiscreteUnif(0,valid_nodes.size()-1);
	treeNode* node_ptr = valid_nodes[index];
	qAnew = qAnew * 1/double(valid_nodes.size());


	// Find the relatives of the selected node
	treeNode* aunt_ptr = NULL;
	treeNode* sib_ptr = NULL;

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

	
	// Find Pr(A|G)
	probA = probP5( node_ptr, aunt_ptr, sib_ptr, theta, rho, first_markers,
					ref_marker, locations, hap_probs, allele_probs );

	
	// Compute the update probabilities Q( A | tilde(A) ). Note that this
	// does not involve simulating new values yet. 
	qA = qA * updateR( 0, aunt_ptr, rho, first_markers, locations,
					   ref_marker, 6, node_ptr, sib_ptr);
	qA = qA * updateR( 0, node_ptr, rho, first_markers, locations,
					   ref_marker, 7, aunt_ptr );
	qA = qA * updateR( 0, sib_ptr->parent, rho, first_markers, locations,
					   ref_marker, 3 );
	qA = qA * updateS( 0, sib_ptr->parent, first_markers, theta, allele_probs,
					   hap_probs );

				
	// Store the old values for the nodes that will change
	treeNode copy_node = *(node_ptr);
	treeNode copy_sib_parent = *(sib_ptr->parent);
	treeNode copy_aunt = *(aunt_ptr);
	treeNode copy_new_parent = *(aunt_ptr->parent);
					
	
	// Change the tree as necessary
	topologyChangeP5(node_ptr, aunt_ptr, sib_ptr);

	
	// Generate new r values for K_i (node_ptr)
	qAnew = qAnew * updateR( 1, node_ptr, rho, first_markers, locations,
							 ref_marker, 6, aunt_ptr, sib_ptr );

	// This first update will modify the z values for K_i's parent, however
	// these values should not change. So I am copying back the original values 
	node_ptr->parent->z_R=copy_new_parent.z_R;
	node_ptr->parent->z_L=copy_new_parent.z_L;

	
	// Generate new r values for , K_q (aunt_ptr) and
	// K_p (node_ptr's old parent or sib_ptr->parent) 
	qAnew = qAnew * updateR( 1, aunt_ptr, rho, first_markers, locations,
							 ref_marker, 7, node_ptr );
	qAnew = qAnew * updateR( 1, sib_ptr->parent,rho, first_markers, locations,
							 ref_marker, 3 );

	
	// Generate a new sequence for K_p (sib_ptr->parent)
	qAnew = qAnew * updateS( 1, sib_ptr->parent, first_markers, theta,
							 allele_probs, hap_probs );


	// For the reverse rearrangement, compute the number of nodes that can
	// be moved in tilde(A)
	valid_nodes.clear();
	findValidNodesP5(valid_nodes, node_ptrs[node_ptrs.size()-1]);
	qA = qA * 1/double(valid_nodes.size());

	
	// Compute the term Pr(tilde(A)|G)
	probAnew = probP5( node_ptr, aunt_ptr, sib_ptr, theta, rho, first_markers,
					   ref_marker, locations, hap_probs, allele_probs );


	// Check that there are no infinite or 0 values as this will cause an
	// undefined MH ratio
	bool errorFlag=0;

	errorFlag = checkValue( qAnew, 1 );
	if (errorFlag==1) {
		/*cout*/ Rcpp::Rcout << "WARNING: Probability of proposed change in minor topology update is infinity" << endl;
	}

	errorFlag = checkValue( qA, 0 );
	if (errorFlag==1) {
		/*cerr*/ Rcpp::Rcout << "ERROR: Probability of proposing current values in minor topology update is 0" << endl;
		throw Rcpp::exception("ERROR: Probability of proposing current values in minor topology update is 0"); //exit(1);
	}

	errorFlag = checkValue( probA, 0 );
	if (errorFlag==1) {
		/*cerr*/ Rcpp::Rcout << "ERROR: Probability of previous data is 0 in minor topology update" << endl;
		throw Rcpp::exception("ERROR: Probability of previous data is 0 in minor topology update"); //exit(1);
	}

	errorFlag = checkValue( probAnew, 1 );
	if (errorFlag==1) {
		/*cout*/ Rcpp::Rcout << "WARNING: Probability of data given proposed changes in minor topology update"
			 << "is infinity. Watch for getting stuck" << endl;
	}



	
	// Compute the MH acceptance probability
	MHratio = min( probAnew / probA * qA / qAnew, double( 1 ) );

	
	// Convert back to the old time configuration since we are not
	// accepting the change
	if ( myunif > MHratio ){
		topologyChangeP5(aunt_ptr, node_ptr, sib_ptr);
		changeR(node_ptr, copy_node);
		changeR(aunt_ptr, copy_aunt);
		changeR(sib_ptr->parent, copy_sib_parent);
		changeS(sib_ptr->parent, copy_sib_parent.seq );
		accept = 0;
	}

	return( accept );
}


// This function is used to find the set of nodes that are eligible for the
// minor topology update. A node is suitable if its parent is older than its
// aunt. 
void findValidNodesP5(vector<treeNode*>& valid_nodes, treeNode* head_node)
{

	if (head_node->child1 != NULL){
		if ((head_node->child1->child1!=NULL)||(head_node->child2->child1!=NULL)){
			
			if (head_node->child1->t > head_node->child2->t){
				valid_nodes.push_back(head_node->child1->child1);
				valid_nodes.push_back(head_node->child1->child2);
			} else {
				valid_nodes.push_back(head_node->child2->child1);
				valid_nodes.push_back(head_node->child2->child2);
			}
			findValidNodesP5(valid_nodes, head_node->child1);
			findValidNodesP5(valid_nodes, head_node->child2);
		}
	}
	
}


// This computes the term proportional to Pr(A|G) for the minor topology update.
//  
double probP5( treeNode* node_ptr, treeNode* aunt_ptr, treeNode* sib_ptr,
			   double theta, double rho, vector<int>& first_markers,
		 	   int ref_marker, vector<double>& locations,
	 		   vector<vector<double> >& hap_probs,
  			   vector<vector<double> >& allele_probs)
{
	double prob = 1;

	// Compute the probability of the sequences and the r's
	treeNode* nodes[] = {node_ptr, aunt_ptr, sib_ptr, node_ptr->parent};

	for (int i=0; i<4; i++){
		prob = prob * probSi( nodes[i], theta, first_markers, allele_probs,
						     hap_probs );
		if ( nodes[i]!= sib_ptr ){
			prob = prob * probRi( nodes[i], rho, first_markers, locations,
								 ref_marker );
		}
	}
	
	return( prob );
	
}


// This function is used to rearrange the nodes for the minor topology
// rearrangement
void topologyChangeP5( treeNode* node_ptr, treeNode* aunt_ptr,
					   treeNode* sib_ptr)
{

	// Change the ptrs of the parents
	if (node_ptr->parent->child1 == node_ptr){
		node_ptr->parent->child1 = aunt_ptr;
	} else {
		node_ptr->parent->child2 = aunt_ptr;
	}

	if (aunt_ptr->parent->child1 == aunt_ptr){
		aunt_ptr->parent->child1 = node_ptr;
	} else {
		aunt_ptr->parent->child2 = node_ptr; 
	}
	
	// Change the ptrs of the children
	aunt_ptr->parent = node_ptr->parent;
	node_ptr->parent = node_ptr->parent->parent;

	// Change the branch lengths
	node_ptr->b = node_ptr->parent->t - node_ptr->t;
	aunt_ptr->b = aunt_ptr->parent->t - aunt_ptr->t;

	// Re-evaluate the z values for the moved node
	aunt_ptr->parent->z_L = max( aunt_ptr->r_L, sib_ptr->r_L );
	aunt_ptr->parent->z_R = max( aunt_ptr->r_R, sib_ptr->r_R );
	
}

