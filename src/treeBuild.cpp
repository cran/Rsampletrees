/*-----------------------------------------------------------------------------

  treeBuild.cpp - This contains the functions for creating a tree and
                  basic manipulation functions on that tree. Note that I
                  expect to replace my tree code with a class version in the 
		  		  future and thus I have created some basic node manipulation
		 		  functions

  Kelly Burkett; March 2007

  Modification Notes:

  ---------------------------------------------------------------------------*/

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <Rcpp.h>

#include "nodefunctions.h" 
#include "readFiles.h" 
#include "rng.h" 
#include "UPGMA.h"
#include "treeBuild.h"
#include "proposal.h"


//***************************************************************************//
// TREE BUILD RELATED FUNCTIONS
//***************************************************************************//

//This function is used by createInitialTree to join internal nodes
treeNode* joinNodes(int id, double time, string seq, vector<int> max_length,
					treeNode* child1, treeNode* child2,
					vector<double> locations, int focal_point, vector<int> first_markers,
					double rho)
{

	treeNode *node_ptr = createNode(id, max_length[0],max_length[1],seq);
	int leftz=0, rightz=0, leftr=0, rightr=0;
	int j=0;

	assignChildren(node_ptr, child1, child2);
	assignParent(child1, node_ptr);
	assignParent(child2, node_ptr);
	changeT(node_ptr, time);


	//Simulate new r values for child 1 (must wait until parent created
	// to get branch time
//	int flag=0;
	updateR(1, child1, rho, first_markers, locations, focal_point, 4);
	updateR(1, child2, rho, first_markers, locations, focal_point, 4);


	//Replace parents parental sequence based on z values

	// Find the locus corresponding to z in the parent
	leftz=first_markers[0]-node_ptr->z_L+1;
	rightz=first_markers[1]+node_ptr->z_R-1;


	if (node_ptr->seq[leftz+1]=='-'){ // Sequence in offspring was shorter on left
		if (child1->r_L > child2->r_L){ // Get sequence from child1

			leftr=first_markers[0]-child2->r_L+1;
			for ( j=leftr; j>leftz; j--){
				node_ptr->seq[j] = child1->seq[j];
			}

		} else { // Get sequence from child2

			leftr=first_markers[0]-child1->r_L+1;
			for ( j=leftr; j>leftz; j--){
				node_ptr->seq[j] = child2->seq[j];
			}

		}
	}

	if (node_ptr->seq[rightz-1]=='-'){ // Sequence in offspring was shorter on right
		if (child1->r_R > child2->r_R){ // Get sequence from child1

			rightr=first_markers[1]+child2->r_R-1;
			for ( j=rightr; j<rightz; j++){
				node_ptr->seq[j] = child1->seq[j];
			}

		} else { // Get sequence from child2

			rightr=first_markers[1]+child1->r_R-1;
			for ( j=rightr; j<rightz; j++){
				node_ptr->seq[j] = child2->seq[j];
			}

		}
	}

	for (j=leftz; j>=0; j--){
		node_ptr->seq[j]='-';
	}
	for (j=rightz; j<(int)node_ptr->seq.size(); j++){
		node_ptr->seq[j]='-';
	}

	return (node_ptr);

}

// This creates the initial tree assuming that no initial tree is provided
// by the user. The internal nodes are connected together based
// on the distance between their respective sequences. The vector of sequence
// information is passed to the function. For each element of this vector, a
// subtree is first made for all the individuals containing the same sequence,
// and a vector of pointers to the heads of these subtrees is kept. Coalescence
// times are then assigned within each subtree to all but the head nodes of the
// subtrees. A UPGMA-like algorithm is then used to determine in which order the
// internal nodes are joined. As these nodes are created, they are assigned a
// coalescence time, id and given the sequence from the subtree with more
// branches at this point.
treeNode* createInitialTree( vector<string> seqs, vector<int> first_markers,
							 vector<double> locations, int focal_point,
							 double rho)
{
	int i=0;

	// Generate the vector of exponential times to be assigned to the internal
	// nodes. T stores the time since the present and is initially set to time
	// of the first coalescent event.
	int num_seq = seqs.size();
    vector<double> times ( seqs.size()-1 );

	for (i=0; i<(int)times.size(); i++){

		// While there are num_seq-i sequences left determine the exponential
		// parameter and generate the times from it.
		times[i] = generateExp((num_seq-i) * (num_seq-i-1) / 2);
	}

	double T=times[0];
	times.erase(times.begin());

	// Set the ID of the first node
	int id = seqs.size()+1;

	// Figure out the maximum values for z and r variables. All tip nodes
	// have z values set to their maximums
	vector<int> max_extent(2);

	max_extent[0] = first_markers[0]+2;
	max_extent[1] = seqs[1].size()-first_markers[1]+1;

	if (first_markers[0] == -1){ // Focal point is before the 1st
		max_extent[0] = 0;
	}
	if (first_markers[1] == (int)seqs[1].size()){// Focal point is after the last
		max_extent[1] = 0;
	}

	// Create a vector of nodeptrs to tip nodes
	vector<treeNode*> node_ptrs(2*seqs.size()-1);

	for (i=0; i<(int)seqs.size(); i++){
		node_ptrs[i] = createNode( i+1,max_extent[0],max_extent[1],seqs[i]);
	}

	// Summarize the sequence information in a struct
	vector<SeqData> seq_info;
	seq_info = uniqueSeq(seqs);


	// To create the tree, we will first be joining all identical sequences
	// together and randomly selecting them based on the frequency of the
	// sequence. To do this, create a vector of sequence frequencies by
	// extracting this information from seq_info. The sequence corresponding
	// to the count will be in seq_info[i].seq_type. Total will be
	// required to keep track of the number of sequences left to coalesce and
	// num_above_1 keeps track of the number of values in seq_freq that are
	// greater than 1. If there is only one sequence of a particular type,
	// have the pointer for that subtree point to that node.
	// For the distance matrix, we will need to know what the unique
	// sequences are (uniqe_seqs) and how many samples have that seq (weights). 
	vector<double> seq_freq(seq_info.size());
	double new_total=0, total = 0;
	int num_above_1 = 0;
	vector<treeNode*> subtree_head_ptrs(seq_info.size());
	vector<string> unique_seqs;
	vector<int> weights;

	for (i=0; i<(int)seq_freq.size(); i++){
		unique_seqs.push_back(seq_info[i].seq_type);
		weights.push_back(seq_info[i].num_of_type);
		if (seq_info[i].num_of_type>1){
			seq_freq[i] = seq_info[i].num_of_type;
			total += seq_freq[i];
			num_above_1++;
			subtree_head_ptrs[i] = NULL;
		} else {
			seq_freq[i]=0;
			subtree_head_ptrs[i] = node_ptrs[seq_info[i].seq_ID[0]-1];
		}
	}

	for (i=0; i<(int)seq_freq.size(); i++){
		seq_freq[i] = seq_freq[i]/total;
	}



	// This is the loop for joining the identical nodes into subtrees
	// with head nodes stored in the vector subtree_head_ptrs
	int current_seq = 0;
	int first=0, second=0, first_ID=0, second_ID=0;
	vector<int>::iterator it_int;

	while (num_above_1>0){

		// Sample a sequence type for joining
		current_seq = samplePMF(seq_freq);

		// Randomly sample two IDs from the set of IDs having the sequence
		// selected.
		first = generateDiscreteUnif(0,seq_info[current_seq].seq_ID.size()-1);
		first_ID = seq_info[current_seq].seq_ID[first];
		it_int = seq_info[current_seq].seq_ID.begin();
		seq_info[current_seq].seq_ID.erase(it_int+first);

		second = generateDiscreteUnif(0,seq_info[current_seq].seq_ID.size()-1);
		second_ID = seq_info[current_seq].seq_ID[second];
		it_int = seq_info[current_seq].seq_ID.begin();
		seq_info[current_seq].seq_ID.erase(it_int+second);

		// Join the two sequences at a new node
		node_ptrs[id-1] = joinNodes( id,T,seq_info[current_seq].seq_type,
									 max_extent,node_ptrs[first_ID-1],node_ptrs[second_ID-1],
									 locations, focal_point, first_markers, rho );
		subtree_head_ptrs[current_seq] = node_ptrs[id-1];


		// Add ID of this node to set of IDs having this sequence
		it_int=seq_info[current_seq].seq_ID.begin();
		seq_info[current_seq].seq_ID.push_back(id);

		// Increase the time and id so these values are ready for the next iteration.
		// Get information for re-weighting
		T += times[0];
		times.erase(times.begin());
		id += 1;

		if (seq_info[current_seq].seq_ID.size()==1){
			num_above_1--;
			seq_freq[current_seq]=0;
			new_total = total-2;
		} else {
			new_total = total-1;
			seq_freq[current_seq] = (seq_freq[current_seq]*total-1)/new_total;
		}
		
		// Re-weight frequency vector
		if (num_above_1 != 0){
			for (i=0; i<(int)seq_freq.size();i++){
				if (i!=current_seq){
					seq_freq[i] = seq_freq[i]*total/new_total;
				}
			}
			total=new_total;
		}

	}

	// Now connect the subtrees together using a UPGMA approach
	
	// First step: Set up distance matrix,
	vector<double> colvecs(seq_info.size());
	vector<vector<double> > distance_mat(seq_info.size(), colvecs);

	computeDistMat(distance_mat, unique_seqs);


	// This is the main loop for connecting the subtrees. Each of the elements
	// in subtree_head_ptrs stores a pointer to the head node of that subtree.
	// From previously, id and T should be set up for the next node to be
	// created
	int pair[2]={0,1};
	vector<treeNode*>::iterator it_treeNode;
	string seq;

	while (distance_mat.size()>=2){
		
		findMinPair(pair, distance_mat);
		
		//cout << "Merge pair: " << pair[0]+1 << "  " << pair[1]+1 << endl;
		//cout << unique_seqs[pair[0]] << " " << unique_seqs[pair[1]] << endl;
		//cout << ",,";

		//for (int l=0; l<distance_mat.size(); l++){
		//	cout  << l+1 << ",";
		//}
		//cout <<endl;
		//for (int l=0; l<distance_mat.size(); l++){
		//	cout << unique_seqs[l] << "," << l+1 << ",";
		//	for (int m=0; m<distance_mat[0].size(); m++){
		//		cout << distance_mat[l][m] << ",";
		//	}
		//	cout << endl;
		//}
		//cout << endl;

		// Set up the sequence information at the head node. Choose the
		// sequence that is more common. Put that sequence information into the
		// pair[0]th element of the sequence vector and delete the information
		// at the pair[1]st element since these two sequences are now merged.
		if ( weights[pair[0]] > weights[pair[1]] ){
			seq = unique_seqs[pair[0]];
		} else {
			seq = unique_seqs[pair[1]];
			unique_seqs[pair[0]]=unique_seqs[pair[1]];
		}
		unique_seqs.erase(unique_seqs.begin()+pair[1]);

		// Join nodes for the minimum pair
		subtree_head_ptrs[pair[0]] = joinNodes( id, T, seq, max_extent,
									        subtree_head_ptrs[pair[0]],
		                                    subtree_head_ptrs[pair[1]],
		                                    locations, focal_point,
		                                    first_markers, rho);

		// Remove row/colum corresponding to minimum pair
		updateDistMat(distance_mat,pair,weights);



		// Update the pointers to the subtrees after joining of minimum pair
		it_treeNode = subtree_head_ptrs.begin();
		subtree_head_ptrs.erase(it_treeNode+pair[1]);

		// Update the number of sequences making up each row/column of
		// th distance matrix
		if (distance_mat.size()!= 1){
			weights[pair[0]] += weights[pair[1]];
			weights.erase(weights.begin()+pair[1]);
			T += times[0];
			times.erase(times.begin());
			id ++;
		}

	}

	// current_head_ptr now points to the MRCA for the whole tree.
	// The r values for this node are 0, since it has no parent
	subtree_head_ptrs[pair[0]]->r_L=0;
	subtree_head_ptrs[pair[0]]->r_R=0;

	return (subtree_head_ptrs[pair[0]]);

}


// This creates the initial tree assuming that no initial tree is provided
// by the user. The tree created is a random one, so that similar sequences
// are not necessarily connected to each other first. 
treeNode* createRandomTree( vector<string> seqs, vector<int> first_markers,
							vector<double> locations, int focal_point,
							double rho)
{
	unsigned int i=0;

	// Generate the vector of exponential times to be assigned to the internal
	// nodes. T stores the time since the present and is initially set to time
	// of the first coalescent event.
	int num_seq = seqs.size();
	vector<double> times ( seqs.size()-1 );

	for (i=0; i<times.size(); i++){

	// While there are num_seq-i sequences left determine the exponential
	// parameter and generate the times from it.
		times[i] = generateExp((num_seq-i) * (num_seq-i-1) / 2);
	}

	double T=times[0];
	times.erase(times.begin());

	// Set the ID of the first node
	int id = seqs.size()+1;

	// Figure out the maximum values for z and r variables. In initial treee
	// will have all r/z values set to the maximum
	vector<int> max_extent(2);

	max_extent[0] = first_markers[0]+2;
	max_extent[1] = seqs[1].size()-first_markers[1]+1;

	if (first_markers[0] == -1){ // Focal point is before the 1st
		max_extent[0] = 0;
	}
	if (first_markers[1] == (int)seqs[1].size()){// Focal point is after the last
		max_extent[1] = 0;
	}

	// Create a vector of nodeptrs to tip nodes
	vector<treeNode*> node_ptrs(seqs.size());

	for (i=0; i<seqs.size(); i++){
		node_ptrs[i] = createNode(i+1,max_extent[0],max_extent[1],seqs[i]);
	}

	// Summarize the sequence information in a struct
	vector<SeqData> seq_info;
	seq_info = uniqueSeq(seqs);

	// Create the tree
	int first=0, second=0;
	
	while (node_ptrs.size()>1){
		
		 first = generateDiscreteUnif(0,node_ptrs.size()-1);
		 do
			 second = generateDiscreteUnif(0,node_ptrs.size()-1);
		 while ( first == second );

		 // Join the first and second node and replace the first with
		 // its parent. Then delete the second.
		 node_ptrs[first]=joinNodes( id, T, node_ptrs[first]->seq, max_extent,
									 node_ptrs[first], node_ptrs[second],
									 locations, focal_point, first_markers, rho);
		 node_ptrs.erase(node_ptrs.begin()+second);
		 
		 // Get ready for the next loop
		 id++;
		 T+=times[0];
		 if (times.size()!=0){
		 	times.erase(times.begin());
		 } else {
			 if (node_ptrs.size()!=1){
				 /*cerr*/ Rcpp::Rcout << "PROBLEM IN RANDOM TREE BUILD";
				 throw Rcpp::exception("PROBLEM IN RANDOM TREE BUILD"); //exit(1);
			 }
		 }

	}
	changeR("left",node_ptrs[0],0);
	changeR("right",node_ptrs[0],0);

	return(node_ptrs[0]);

}



//***************************************************************************//
// TREE BUILD FUNCTIONS WHEN GIVEN THE INITIAL TREE IN FILE 
//***************************************************************************//


// Generate the initial tree from information in the file given. The tree
// should be given in a tabular format in the file and is likely what was
// output by a previous run of this function. The functions 
// readInitialTree, fillTree and fillNode are used to do this.
treeNode* readInitialTree ( ifstream& initial_file )
{
	treeNode *head_ptr = createNode();

	fillTree( head_ptr, NULL, initial_file );

	return ( head_ptr );
}

// Recursive function to traverse the tree and fill in nodes as we go
void fillTree ( treeNode* tree_ptr, treeNode* parent_ptr, ifstream& file )
{

	fillNode ( tree_ptr, parent_ptr, file );

	if ( tree_ptr->child1 != NULL ){
		fillTree ( tree_ptr->child1, tree_ptr, file );
	}
	if ( tree_ptr->child2 != NULL ){
		fillTree ( tree_ptr->child2, tree_ptr, file );
	}

}

// Function to actually fill in the node based on the information from the
// input file.
void fillNode ( treeNode* tree_ptr, treeNode* parent_ptr, ifstream& file )
{
	int parent_id, child1_id;

	file >> tree_ptr->id;
	file >> parent_id;
	tree_ptr->parent = parent_ptr;
	file >> child1_id;
	if ( child1_id == 0 ){
		tree_ptr->child1 = NULL;
		tree_ptr->child2 = NULL;
	} else { // they will be filled in later
		tree_ptr->child1 = createNode();
		tree_ptr->child2 = createNode();
	}
	file >> child1_id;
	file >> tree_ptr->t;
	file >> tree_ptr->b;
	file >> tree_ptr->z_L;
	file >> tree_ptr->r_L;
	file >> tree_ptr->z_R;
	file >> tree_ptr->r_R;
	file >> tree_ptr->seq;

}


//***************************************************************************//
// TREE REMOVAL RELATED FUNCTIONS
//***************************************************************************//


// These two functions are used to delete the tree. Since there are two ways
// to access the tree, either by node_ptrs or with the head_node, there are
// two ways to delete it. 
void deleteTree ( vector<treeNode*>& node_ptrs )
{
	unsigned int i=0;
	for (i=0; i<node_ptrs.size(); i++){
		delete node_ptrs[i];
		node_ptrs[i] = NULL;
	}
}



void deleteTree ( treeNode* node_ptr )
{
	if (node_ptr == NULL){
		return;
	}
	deleteTree(node_ptr->child1);
	deleteTree(node_ptr->child2);
	if (node_ptr->parent !=NULL ){
		if (node_ptr->parent->child1 == node_ptr){
			node_ptr->parent->child1 = NULL;
		} else {
			node_ptr->parent->child2 = NULL;
		}
	}
	//cout << "Deleting node with ID" << node_ptr->id << endl;
	delete node_ptr;
	node_ptr = NULL;
}


//***************************************************************************//
// TREE COPY FUNCTION
//***************************************************************************//


// This is the function that copies a tree that is pointed to by node_ptr.
// It traverses the tree until it reaches the bottom of the tree (on the
// left hand side) and creates copies of the nodes as it traverses back
// up the tree. This function sets the value of one element of each node,
// the pointer to the node's parent. This is done this way since the
// parent node cannot be set by the copyNode function.
treeNode* copyTree(treeNode* node_ptr)
{
	treeNode *left_ptr, *right_ptr, *current_ptr;

	if (node_ptr == NULL){
		return NULL;
	} else {
		left_ptr = copyTree(node_ptr -> child1);
		right_ptr = copyTree(node_ptr ->child2);
		current_ptr = copyNode(node_ptr, left_ptr, right_ptr);
		if ( current_ptr -> child1 != NULL ){
			current_ptr->child1->parent = current_ptr;
			current_ptr->child2->parent = current_ptr;
		}
		return(current_ptr);
	}
}


//***************************************************************************//
// PRINT RELATED FUNCTIONS
//***************************************************************************//

// Wrapper function to print the tree. It calls one of the printTreeStruc
// functions, depending on whether the times are required or not.
void printTreeStruc(treeNode* node_ptr,int iteration, bool printTimes, ostream& outfile){
	
	if (printTimes == 0){
		outfile << iteration << " ";
		printTree(node_ptr, outfile);
		outfile << ";" << endl;
	}
	if (printTimes == 1){
		outfile << iteration << " ";
		printTree(node_ptr, printTimes, outfile);
		outfile << ";" << endl;
	}
}


// This function is used to print out the tree information in NEWICK format
// minus the branch lengths. It prints out information about the topology
// of the tree only.
void printTree(treeNode* node_ptr, ostream& outfile){

	outfile  << "(";

	if ((node_ptr->child1->child1 == NULL)&&(node_ptr->child2->child1 == NULL)){

		// Left and Right grandchildren of current node are NULL (at bottom)
		outfile << node_ptr->child1->id << "," << node_ptr->child2->id;

	} else if ((node_ptr->child1->child1 != NULL)&&(node_ptr->child2->child1 == NULL)){

		// Only Right grandchild is NULL. So at bottom on right side of tree, but
		// not left, so print the left side sub-tree
		printTree(node_ptr->child1, outfile);
		outfile << "," << node_ptr->child2->id;

	} else if ((node_ptr->child1->child1 == NULL)&&(node_ptr->child2->child1 != NULL)){

		// Only Left grandchild is NULL. So at bottom on left side of tree, but
		// not right, so print the right side sub-tree
		outfile << node_ptr->child1->id << ",";
		printTree(node_ptr->child2, outfile);

	} else {

		// Subtrees on both left and right hand side to print
		printTree(node_ptr->child1, outfile);
		outfile << ",";
		printTree(node_ptr->child2, outfile);

	}

	outfile << ")";

}


// This is the same as the above, but with the option of printing times out as
// well. The times can be used by PHYLIP - drawgram to draw a phylogenetic
// tree with branches equal to the branch length "b"
void printTree(treeNode* node_ptr, bool printTimes, ostream& outfile){

	outfile << "(";

	if ((node_ptr->child1->child1 == NULL)&&(node_ptr->child2->child1 == NULL)){

		// Left and Right grandchildren of current node are NULL (at bottom)
		outfile << node_ptr->child1->id << ":" << node_ptr->child1->b;
		outfile <<	",";
		outfile << node_ptr->child2->id << ":" << node_ptr->child2->b;

	} else if ((node_ptr->child1->child1 != NULL)&&(node_ptr->child2->child1 == NULL)){

		// Only Right grandchild is NULL. So at bottom on right side of tree, but
		// not left, so print the left side sub-tree
		printTree(node_ptr->child1, 1, outfile);
		outfile << ":" << node_ptr->child1->b;
		outfile << "," << node_ptr->child2->id;
		outfile << ":" << node_ptr->child2->b;


	} else if ((node_ptr->child1->child1 == NULL)&&(node_ptr->child2->child1 != NULL)){

		// Only Left grandchild is NULL. So at bottom on left side of tree, but
		// not right, so print the right side sub-tree
		outfile << node_ptr->child1->id << ":" << node_ptr->child1->b << ",";
		printTree(node_ptr->child2, 1, outfile);
		outfile << ":" << node_ptr->child2->b;

	} else {

		// Subtrees on both left and right hand side to print
		printTree(node_ptr->child1,1, outfile);
		outfile << ":" << node_ptr->child1->b;
		outfile << ",";
		printTree(node_ptr->child2,1, outfile);
		outfile << ":" << node_ptr->child2->b;

	}

	outfile << ")";

}


// This function is used to print out all the relevant information at each of
// the nodes. One line is given per node. The order is:
// ID Parent child1 child2 t b z_L r_L seq

void printTreeFull(treeNode *node_ptr, ostream& outfile)
{
		if ( node_ptr != NULL ) {  // (Otherwise, there's nothing to print.)

			outfile << node_ptr->id << " ";

			// We are at the top of the tree so there is no parent
			if (node_ptr->parent == NULL){
				outfile << "0 ";
			} else {
				outfile << node_ptr->parent->id << " ";
			}

			// We are at the bottom of the tree so there are no children
			if (node_ptr->child1 == NULL){
				outfile << "0 0 ";
			} else {
				outfile << node_ptr->child1->id << " " << node_ptr->child2->id << " ";
			}

			// Print out all info required
			outfile << node_ptr->t << " " << node_ptr->b << " " ;
			outfile << node_ptr->z_L << " " << node_ptr->r_L << " ";
			outfile << node_ptr->z_R << " " << node_ptr->r_R << " ";
			outfile << node_ptr->seq << endl;

			// Print the left and right subtree respectively
			printTreeFull( node_ptr->child1 , outfile);
			printTreeFull( node_ptr->child2 , outfile);
		}
}





//***************************************************************************//
// MISCELLANEOUS FUNCTIONS
//***************************************************************************//




// This function traverses the tree to set pointers to each of the nodes in
// the tree. These are stored in the vector node_ptrs. This can be used
// for sampling nodes if sampling is done uniformly.
void getNodePtrs(vector<treeNode*>& node_ptrs, treeNode *root_ptr)
{
	int index;

	if ( root_ptr != NULL ) {
		getNodePtrs(node_ptrs, root_ptr->child1);
		getNodePtrs(node_ptrs, root_ptr->child2);
		index = root_ptr->id-1;
		node_ptrs[index] = root_ptr;
	}
}


// This function takes the vector of node_ptrs and orders it based on
// increasing time of internal nodes. The first n nodes do not change order
// as they correspond to the tip nodes at present. The second version of 
// this function can be used if a known subvector of node_ptrs has changed.
// The indices definining the subvector are passed and only those 
// elements between the two are sorted.
void sortNodePtrs(vector<treeNode*>& node_ptrs)
{
	vector<treeNode*>::iterator it;
	treeNode* copy_node = NULL;
	bool found=0;
	int i=0, j=0, k=0;
	int n = (node_ptrs.size()+1)/2;

	i=n;
	while (i < (int)node_ptrs.size()){

		// If the time of the current node is less than the time at the
		// previous node then a re-ordering is needed.
		if ( node_ptrs[i]->t < node_ptrs[i-1]->t ){

			// Make a copy of the node to be moved and delete it from the
			// vector
			copy_node = node_ptrs[i];
			it = node_ptrs.begin()+i;
			node_ptrs.erase(it);

			// Now find where the node must be inserted. Since t_i < t_{i-1} we
			// must decrease j until we find a node with t_j <= t_i
			// and t_{j+1}>t_i. Then the node is inserted at j+1 and all elements
			// between the old location and the new location get a new id.
			j=i-1;
			while (found==0){
				if (node_ptrs[j]->t <= copy_node->t){
					found = 1;
					it = node_ptrs.begin()+j+1;
					node_ptrs.insert(it,copy_node);
				} else {
					j--;
				}
			}
			for (k=j; k<i+1; k++){
				node_ptrs[k]->id = k+1;
			}
			found=0;
		}

		// Ensure that the id at this node is correct
		node_ptrs[i]->id = i+1;
		i++;

	}

	for (i=0; i<(int)node_ptrs.size()-1; i++){
		if (node_ptrs[i]->t > node_ptrs[i+1]->t){
			/*cout*/ Rcpp::Rcout << "WARNING: sort failure of times for id " << i << endl;
		}
		if (node_ptrs[i]->id > node_ptrs[i+1]->id){
			/*cout*/ Rcpp::Rcout << "WARNING: sort failure of id's for id " << i << endl;
		}
	}

}


void sortNodePtrs(vector<treeNode*>& node_ptrs, int min_index, int max_index)
{
	vector<treeNode*>::iterator it;
	treeNode* copy_node = NULL;
	bool found=0;
	int i=0, j=0, k=0;

	i=min_index;
	
	while (i < max_index ){

			// If the time of the current node is less than the time at the
			// previous node then a re-ordering is needed.
			if ( node_ptrs[i]->t < node_ptrs[i-1]->t ){


				// Make a copy of the node to be moved and delete it from the
				// vector
				copy_node = node_ptrs[i];
				it = node_ptrs.begin()+i;
				node_ptrs.erase(it);

				// Now find where the node must be inserted. Since t_i < t_{i-1} we
				// must decrease j until we find a node with t_j <= t_i
				// and t_{j+1}>t_i. Then the node is inserted at j+1 and all elements
				// between the old location and the new location get a new id.
				j=i-1;
				while (found==0){
					if (node_ptrs[j]->t <= copy_node->t){
						found = 1;
						it = node_ptrs.begin()+j+1;
						node_ptrs.insert(it,copy_node);
					} else {
						j--;
					}
				}
				// Ensure that all previous nodes ids are correct
				for (k=j; k<i+1; k++){ 
					node_ptrs[k]->id = k+1;
				}
				found=0;
			}

			// Ensure that the id at this node is correct
			node_ptrs[i]->id = i+1;
			i++;

		}


	for (i=0; i<(int)node_ptrs.size()-1; i++){
		if (node_ptrs[i]->t > node_ptrs[i+1]->t){
			/*cout*/ Rcpp::Rcout << "WARNING: sort failure of times for id " << i << endl;
		}
		if (node_ptrs[i]->id > node_ptrs[i+1]->id){
			/*cout*/ Rcpp::Rcout << "WARNING: sort failure of id's for id " << i << endl;
		}

	}
	
}






