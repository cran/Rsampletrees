/*-----------------------------------------------------------------------------

  nodefunctions.h - This is the header file for nodefunctions.cpp, which 
		            contains the functions to create new nodes and do basic
		            manipulation of the values of the node. I have set it
		            up so it might be possible to have this be a class in
		            the future.

  Kelly Burkett; July 2009

  ---------------------------------------------------------------------------*/

#ifndef NODE_H
#define NODE_H


#include <string>
#include <vector>

using namespace std;

// Definition of the node structure for the tree
struct treeNode
{
	int z_L, z_R, r_L, r_R; 	
	int id;
	double b, t;
	string seq;
	treeNode *child1, *child2, *parent;
};


// These two functions are used to create a new node. The first initializes
// the values to null values, the second can be used when the node is to
// be initialized to known values for id,z,and seq
treeNode* createNode();
treeNode* createNode(int id, int zL, int zR, string seq);

// These two functions set the child/parent pointers for a node
void assignParent(treeNode* node_ptr, treeNode* parent_ptr);
void assignChildren(treeNode* node_ptr, treeNode* child_one, treeNode* child_two);

// This function determines the sib of the node that is passed
treeNode* findsib(treeNode* node_ptr);

// The following functions are used to modify the values of the treeNode
// elements. 
void changeT(treeNode* node_ptr, double new_t);
void changeR(string side, treeNode* node_ptr, int new_r);
void changeR(treeNode* node_ptr, treeNode* copy_ptr);
void changeR(treeNode* node_ptr, treeNode copy_node);
void changeSij( treeNode* node_ptr, char new_sij, int index );
void changeS( treeNode* node_ptr, string new_seq );

// The following two functions make a copy of a node and may include new values
// for the child pointers as this function might be called when making a copy
// of a full tree
treeNode* copyNode(treeNode* node_ptr, treeNode* left_ptr, treeNode* right_ptr);
treeNode* copyNode(treeNode* node_ptr);

#endif

