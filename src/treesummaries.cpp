/*-----------------------------------------------------------------------------

  treesummaries.cpp - This file contains the functions that summarize 
		      		  information about the tree. These functions are not
		      		  currently being used as these summaries are being done
		      		  in R

  Kelly Burkett; July 2009

  Modification Notes:

  ---------------------------------------------------------------------------*/

#include <vector>
#include <algorithm>

#include "treesummaries.h"
#include "nodefunctions.h"


// A quick function that can be used to count nodes in a binary tree.
// I got it off of http://math.hws.edu/eck/cs225/s03/binary_trees/
int countNodes( treeNode *root_ptr )
{
           // Count the nodes in the binary tree to which
           // root points, and return the answer.
	if ( root_ptr == NULL )
		return 0;  // The tree is empty.  It contains no nodes.
	else {
		int count = 1;   // Start by counting the root.
		count += countNodes(root_ptr->child1);  // Add the number of nodes
                                            //     in the left subtree.
		count += countNodes(root_ptr->child2); // Add the number of nodes
                                            //    in the right subtree.
		return count;  // Return the total.
	}
}



// This function returns the level of the MRCA of two nodes.
int get_mrca_level(treeNode* node_ptr_1, treeNode* node_ptr_2, int level)
{

	if (node_ptr_1->parent == node_ptr_2->parent){
		return(level);
	} else{
		if (node_ptr_1->parent->t < node_ptr_2->parent->t){
			level = get_mrca_level(node_ptr_1->parent,node_ptr_2, level+1);
		} else {
			level = get_mrca_level(node_ptr_1,node_ptr_2->parent, level+1);
		}
		return(level);
	}
}


// This function returns the time since the MRCA of particular pair
// of nodes
double get_mrca_time(const treeNode* node_ptr_1, const treeNode* node_ptr_2)
{
	double time=0;

	if (node_ptr_1->parent != node_ptr_2->parent)
	{
		if (node_ptr_1->parent->t < node_ptr_2->parent->t){
			time = get_mrca_time(node_ptr_1->parent,node_ptr_2);
		} else {
			time = get_mrca_time(node_ptr_1,node_ptr_2->parent);
		}
		return(time);
	} else {
		return(node_ptr_1->parent->t);
	}
}


// This function returns the average time to MRCA of all pairs of
// tips
double average_mrca_time(const vector<treeNode*>& node_ptrs)
{
	int n = (node_ptrs.size()+1)/2;
	int i=0,j=0;
	double time=0, alltime=0;

	for (i=0; i<n; i++){
		for (j=0; j<i; j++){
			time = get_mrca_time(node_ptrs[i],node_ptrs[j]);
			alltime = alltime + time;
		}
	}
	
	return ( 2*alltime/(n*(n-1)) );
}




// This function fills in a symmetric matrix with the time to mrca for pairs
// of sequences indexed by the row and column. Diagonal is 0.
void average_time_mat( const vector<treeNode*>& node_ptrs,
		       vector<vector<double> >& time_mat)
{
	int n = (node_ptrs.size()+1)/2;
	int i=0,j=0;
	double time=0;
	
	for (i=0; i<n; i++){
		for (j=0; j<i; j++){
			time = get_mrca_time(node_ptrs[i],node_ptrs[j]);
			time_mat[i][j] = time;
			time_mat[j][i] = time;
		}
	}
}



/* --------------------------------------------------------------------------
   The following sets of functions all relate to getting partitions of the 
   tree.
   -------------------------------------------------------------------------*/


// This function returns the partitions of a given tree
vector<int> getPartitions(vector<vector<int> >& partitions, treeNode* node_ptr)
{
	vector<int> current_vec(0);
	vector<int> second_vec(0);

	if ((node_ptr->child1->child1 == NULL)&&(node_ptr->child2->child1 == NULL)){

		// Left and Right grandchildren of current node are NULL (at bottom)
		current_vec.push_back(node_ptr->child1->id);
		current_vec.push_back(node_ptr->child2->id);

	} else if ((node_ptr->child1->child1 != NULL)&&(node_ptr->child2->child1 == NULL)){

		// Only Right grandchild is NULL. So at bottom on right side of tree, but
		// not left, so print the left side sub-tree
		current_vec = getPartitions(partitions, node_ptr->child1);
		current_vec.push_back(node_ptr->child2->id);

	} else if ((node_ptr->child1->child1 == NULL)&&(node_ptr->child2->child1 != NULL)){

		// Only Left grandchild is NULL. So at bottom on left side of tree, but
		// not right, so print the right side sub-tree
		current_vec = getPartitions(partitions, node_ptr->child2);
		current_vec.push_back(node_ptr->child1->id);

	} else {

		// Subtrees on both left and right hand side to print
		current_vec = getPartitions(partitions, node_ptr->child1);
		second_vec =getPartitions(partitions, node_ptr->child2);

		for (unsigned i=0; i<second_vec.size(); i++){
			current_vec.push_back(second_vec[i]);
		}

	}
	if (node_ptr->parent!=NULL){
		partitions.push_back(current_vec);
	}
	return(current_vec);
}


bool comparePartition(vector<int> partition1, vector<int> partition2){

	unsigned i/*, j*/;

	if (partition1.size() != partition2.size()){
		return(0);
	} else {
		for (i=0; i<partition1.size(); i++){
			if (partition1[i] != partition2[i]){
				return(0);
			}
		}
		return(1);
	}
}


// This functions takes two sets of partitions and finds the distance between
// them. The partition sets must be in the appropriate form
int partitionDist(vector<vector<int> >& set_1, vector<vector<int> > set_2)
{
	unsigned i, j;
	bool isEqual = 0;
	int dist = set_1.size() + set_2.size();
	vector<vector<int> >::iterator it;

	for (i=0; i<set_1.size(); i++){
		for (j=0; j<set_2.size(); j++){
			isEqual = comparePartition(set_1[i], set_2[j]);
			if (isEqual==1){
				dist-=2;
				it = set_2.begin()+j;
				set_2.erase(it);
			}

		}
	}

	return(dist);

}

int partitionDist2(vector<vector<int> > set_1, vector<vector<int> > set_2)
{
	unsigned i=0, j=0;
	bool isEqual = 0;
	int dist = set_1.size() + set_2.size();
	vector<vector<int> >::iterator it;


	while ( i<set_1.size()){
		while ( (j<set_2.size()) && (!isEqual) ){
			isEqual = comparePartition(set_1[i], set_2[j]);
			if (isEqual){
				dist-=2;
				it = set_2.begin()+j;
				set_2.erase(it);
			} else {
				j++;
			}
		}
		if (isEqual){
			i++;
			isEqual = 0;
		} else {
			i++;
		}
		j=0;
	}


	return(dist);

}


// This function is used to convert the partitions into a standardized format:
// - Sorted in ascending order
// - The shortest number of elements are chosen. Eg: with 8 terminal nodes, 1,3,4
//   represents the same partition as 2,5,6,7,8. The former will be the
//   standard way of representing it.
// - For even number of terminal nodes, partitions of equal size exist. In this
//   case, the standard partition will be the one that includes 1. Eg: Of the
//   two partitions 1,4,5,7 and 2,3,6,8, the former will be used.
void standardizePartition(vector<vector<int> >& partitions, int num_seqs)
{
	unsigned index, j=0;
  int i=0;
	vector<int> alternate(0);
	vector<vector<int> >::iterator it;

	for (index=0; index<partitions.size(); index++){

		// Sort the current partition
		sort(partitions[index].begin(), partitions[index].end());
		j=0;
		alternate.clear();

		if (partitions[index].size()>double(num_seqs)/2){
			// The current partition is being represented by its longer
			// equivalent. Change to the shorter


			for (i=1; i<=num_seqs; i++){
				if ((j<partitions[index].size())&&(partitions[index][j]==i)){
					j++;
				} else {
					alternate.push_back(i);
				}
			}

			// Erase the current partition and replace it with the shorter
			// version in the vector of all partitions
			it = partitions.begin()+index;
			partitions.erase(it);
			partitions.insert(it,alternate);

		} else if ( (partitions[index].size()==double(num_seqs)/2) &&
						(partitions[index][0]!=1)){
			// The current partition is the same size as its alternate, but
			// doesn't include 1, so change to the alternate

			for (i=1; i<=num_seqs; i++){
				if (partitions[index][j]==i){
					j++;
				} else {
					alternate.push_back(i);
				}
			}

			// Erase the current partition and replace it with the shorter
			// version in the vector of all partitions
			it = partitions.begin()+index;
			partitions.erase(it);
			partitions.insert(it,alternate);

		}
	}

}


