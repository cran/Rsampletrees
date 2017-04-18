/*----------------------------------------------------------------------------

  treeBuild.h - This is the header file for treeBuild.cpp. The node strucutre
		of the tree is defined in this function
		
  Kelly Burkett; March 2007

  --------------------------------------------------------------------------*/

#ifndef TREEBUILD_H
#define TREEBUILD_H

#include <sstream>
#include <string>
#include <vector>

#include "nodefunctions.h"
//#include "sampletrees.h"

using namespace std;

/*
struct treeOut
{
 vector<int> iteration;
 vector<string> tree;
};
*/


// This function is used to join two internal nodes when creating a tree
treeNode* joinNodes( int id, double time, string seq, vector<int> max_length,
		             treeNode* child1, treeNode* child2, vector<int> locations,
		             int focal_point, vector<double> first_markers,
		             double rho);


// This is the main function to create the full tree given the present day
// sequences. It creates a number of subtrees and then uses the UPGMA method
// to connect the subtrees together by most similar. It returns a pointer
// to the head node.
treeNode* createInitialTree( vector<string> seqs, vector<int> first_markers,
							 vector<double> locations, int focal_point,
							 double rho );
treeNode* createRandomTree( vector<string> seqs, vector<int> first_markers,
							vector<double> locations, int focal_point,
							double rho);


// These functions are used to initialize the tree with information from a file
// that is providid by the user
treeNode* readInitialTree ( ifstream& initial_file );
void fillTree ( treeNode* tree_ptr, treeNode* parent_ptr, ifstream& file );
void fillNode ( treeNode* tree_ptr, treeNode* parent_ptr, ifstream& file );


// These functions delete the created tree, either by deleting the nodes 
// pointed to by node_ptrs or by traversing the tree and deleting bottom-up
void deleteTree ( vector<treeNode*>& node_ptrs );
void deleteTree ( treeNode* node_ptr );


// This function makes a copy of the tree having head node node_ptr
treeNode* copyTree(treeNode* node_ptr);


// This is the function to print out the tree topology in the NEWICK
// format. Will have to add coalescent times in the future. 
void printTree(treeNode* node_ptr, ostream& outfile);
void printTree(treeNode* node_ptr, bool printTimes, ostream& outfile);
void printTreeStruc(treeNode* node_ptr, int i, bool printTimes, ostream& outfile);
void printTreeFull(treeNode *node_ptr, ostream& outfile);
				 
				 
// This fills in node_ptrs, which is a vector of pointers to each of the
// nodes in the tree. I think it will be useful when we need to sample
// from the nodes uniformly for making alterations on the tree. 
void getNodePtrs(vector<treeNode*>& node_ptrs, treeNode *root_ptr);


// These functions sort the vector of node_ptrs so that they are in
// order of increasing age. The updates can change the order of coalescence
// events and thus the need to sort. This is an overloaded function so
// that it only sorts between nodes that are known to have changed in the
// second version. This means that all times must be between that of node at
// min_index and max_index!
void sortNodePtrs(vector<treeNode*>& node_ptrs);
void sortNodePtrs(vector<treeNode*>& node_ptrs, int min_index, int max_index);





#endif 

