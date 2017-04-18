/*-----------------------------------------------------------------------------

  treesummaries.h - This is the header file for treesummaries.cpp,

  Kelly Burkett; July 2009

  ---------------------------------------------------------------------------*/

#include <vector>

#include "nodefunctions.h"


#ifndef SUMMARY_H
#define SUMMARY_H

// Compute the number of nodes that the tree has
int countNodes( treeNode *root_ptr );


// Find the number of generations back until the MRCA of the two nodes 
// pointed to by node_ptr_1 and node_ptr_2
int get_mrca_level( treeNode* node_ptr_1, treeNode* node_ptr_2, int level );

// Find the time back until the MRCA of the two nodes pointed to by
// node_ptr_1 and node_ptr_2
double get_mrca_time( const treeNode* node_ptr_1, const treeNode* node_ptr_2 );


// Find the overall average of the times to MRCA for all pairs of sequences
// observed at present
double average_mrca_time( const vector<treeNode*>& node_ptrs );


// Compute a matrix with the MRCA of the two nodes indexed by the row and 
// column. This is a symmetric matrix and the diagonal is 0
void MRCA_time_mat( const vector<treeNode*>& node_ptrs,
		    vector<vector<double> >& time_mat);


// Functions for the Partition Distance summary measure for the tree topology
vector<int> getPartitions(vector<vector<int> >& partitions, treeNode* node_ptr);
int partitionDist(vector<vector<int> >& set_1, vector<vector<int> > set_2);
bool comparePartition(vector<int> partition1, vector<int> partition2, int num_seqs);
void standardizePartition(vector<vector<int> >& partitions, int num_seqs);
int partitionDist2(vector<vector<int> > set_1, vector<vector<int> > set_2);

#endif

