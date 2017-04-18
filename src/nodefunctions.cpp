/*-----------------------------------------------------------------------------

  nodefunctions.cpp - This file contains the functions relating to the node
                      definitions and functions. I have set this up so that
		              it may be possible to set it up as a class in the
		              future

  Kelly Burkett; July 2009

  Modification Notes:

  ---------------------------------------------------------------------------*/

#include <string>
#include <iostream> //For cerr/cout
#include <cstdlib>
#include <Rcpp.h>

#include "nodefunctions.h"
#include "rng.h"



// This function creates a new node and initializes its values. A pointer to
// this node is returned.
treeNode* createNode()
{
	treeNode *node_ptr = new treeNode;

	if (node_ptr == NULL){
		Rcpp::Rcout << "Error: Insufficient memory." << endl;	//cerr << "Error: Insufficient memory." << endl;
		throw Rcpp::exception("Error: Insufficient memory.");	//exit(1);
		
	}

	node_ptr->z_L = 0;
	node_ptr->z_R = 0;
	node_ptr->r_L = 0;
	node_ptr->r_R = 0;
	node_ptr->id = 0;
	node_ptr->b = 0;
	node_ptr->t = 0;
	node_ptr->seq = "";
	node_ptr->child1 = NULL;
	node_ptr->child2 = NULL;
	node_ptr->parent = NULL;

	return (node_ptr);
}


// This function creates a new node and initializes its values, using the
// values provided where applicable. A pointer to this node is returned.
treeNode* createNode(int id, int zL, int zR, string seq)
{
	treeNode *node_ptr = new treeNode;

	if (node_ptr == NULL){
		Rcpp::Rcout << "Error: Insufficient memory." << endl; //cerr << "Error: Insufficient memory." << endl;
		throw Rcpp::exception("Error: Insufficient memory."); //exit(1);
	}

	node_ptr->z_L=zL;
	node_ptr->z_R=zR;
	node_ptr->r_L=zL;
	node_ptr->r_R=zR;

	node_ptr->id = id;
	node_ptr->b = 0;
	node_ptr->t = 0;
	node_ptr->seq = seq;
	node_ptr->child1 = NULL;
	node_ptr->child2 = NULL;
	node_ptr->parent = NULL;

	return (node_ptr);
}


// The following two functions assign the relationships. The functions
// are not necessarily needed but this may make it easier to convert to
// a class for the treeNode later
void assignParent(treeNode* node_ptr, treeNode* parent_ptr)
{
	node_ptr-> parent = parent_ptr;
}

void assignChildren(treeNode* node_ptr, treeNode* child_one, treeNode* child_two)
{
	node_ptr->child1 = child_one;
	node_ptr->child2 = child_two;

}

// This function finds the sib of the node that is passed
treeNode* findsib(treeNode* node_ptr){

	if (node_ptr->parent->child1==node_ptr){
		return(node_ptr->parent->child2);
	} else {
		return(node_ptr->parent->child1);
	}
}

// This function changes the time since present at a node pointed to by
// node_ptr. When t changes, so do the branch lengths for the parents
// and children of the node.
void changeT(treeNode* node_ptr, double new_t)
{
	node_ptr->t = new_t;

	if (node_ptr->parent != NULL){
		node_ptr->b = node_ptr->parent->t - node_ptr->t;
	}
	node_ptr->child1->b = node_ptr->t - node_ptr->child1->t;
	node_ptr->child2->b = node_ptr->t - node_ptr->child2->t;
}


// I have overloaded the changeR function to handle the case of passing
// a node with old values or a pointer to a node with old values, and
// with changing to a new value that is passed to the function.
void changeR(string side, treeNode* node_ptr, int new_r)
{
	if (side == "left"){
		node_ptr->r_L = new_r;
		if (node_ptr->parent!=NULL){
			node_ptr->parent->z_L = max( node_ptr->parent->child1->r_L,
					     node_ptr->parent->child2->r_L);
		}
	} else {
		node_ptr->r_R = new_r;
		if (node_ptr->parent!=NULL){
		node_ptr->parent->z_R = max( node_ptr->parent->child1->r_R,
					     node_ptr->parent->child2->r_R);
		}
	}
}

void changeR(treeNode* node_ptr, treeNode* copy_ptr)
{
	node_ptr->r_L = copy_ptr->r_L;
	node_ptr->r_R = copy_ptr->r_R;
	node_ptr->parent->z_L = max( node_ptr->parent->child1->r_L,
				     node_ptr->parent->child2->r_L);
	node_ptr->parent->z_R = max( node_ptr->parent->child1->r_R,
				     node_ptr->parent->child2->r_R);
}

void changeR(treeNode* node_ptr, treeNode copy_node)
{
	node_ptr->r_L = copy_node.r_L;
	node_ptr->r_R = copy_node.r_R;
	node_ptr->parent->z_L = max( node_ptr->parent->child1->r_L,
				     node_ptr->parent->child2->r_L);
	node_ptr->parent->z_R = max( node_ptr->parent->child1->r_R,
				     node_ptr->parent->child2->r_R);
}




// These two functions change the sequence information. In the first, 
// an element of the sequence is changed. In the second the whole sequence
// is changed.
void changeSij( treeNode* node_ptr, char new_sij, int index )
{
	node_ptr->seq[index] = new_sij;
}

void changeS( treeNode* node_ptr, string new_seq )
{
	node_ptr->seq = new_seq;
}


// This function is used to copy the contents of one node into a new node.
// However, we don't want to copy the locations of the pointers to the children
// as this would point to the wrong tree. So the locations of the children are
// passed to the function. Note that the parent node is not passed. This is
// because the parent may not yet exist in the way that the function is
// currently being called.
treeNode* copyNode(treeNode* node_ptr, treeNode* left_ptr, treeNode* right_ptr)
{
	treeNode *copy_ptr = createNode();

	if (copy_ptr == NULL){
		Rcpp::Rcout << "Error: Insufficient memory." << endl; // cerr << "Error: Insufficient memory." << endl;
		throw Rcpp::exception("Error: Insufficient memory."); //exit(1);
	}

	*(copy_ptr) = *(node_ptr);
	copy_ptr->child1 = left_ptr;
	copy_ptr->child2 = right_ptr;
	copy_ptr->parent = NULL;

	return( copy_ptr );
}


// In this version the parent and children pointers are given as NULL
treeNode* copyNode(treeNode* node_ptr)
{
	treeNode *copy_ptr = createNode();

	if (copy_ptr == NULL){
		Rcpp::Rcout << "Error: Insufficient memory." << endl; // cerr << "Error: Insufficient memory." << endl;
		throw Rcpp::exception("Error: Insufficient memory."); //exit(1);
	}

	*(copy_ptr) = *(node_ptr);
	copy_ptr->parent = NULL;
	copy_ptr->child1 = NULL;
	copy_ptr->child2 = NULL;

	return( copy_ptr );
}
