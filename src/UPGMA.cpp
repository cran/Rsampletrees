/*----------------------------------------------------------------------------

  UPGMA.cpp - This contains the functions for creating a phylogenetic tree
              using the UPGMA method. See "Inferring Phylogenies" by 
              Felsenstein for a reference. Note that there isn't a function
              here which will create the full tree - that is done in
              treeBuild.cpp - but the functions in here are used. The
              reason is that unlike a typical use of this method, I need to
              initialize values at the internal nodes. It seems easier to
              me to do this at the same time as the tree itself is created.

  Kelly Burkett; March 2007

  Modification Notes:

  --------------------------------------------------------------------------*/

#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <random>

#include "UPGMA.h"

using namespace std; 


// This function returns the filled in distance matrix (distmat). Each element
// of distmat is the distance between the two sequences being compared. The 
// rows and columns represent sequences 1 through seqslen. The distance 
// is defined to be the number of alleles that differ between the
// two sequences.
void computeDistMat(vector<vector<double> >& distmat, const vector<string> seqs){
	
	unsigned int i=0,j=0;
  int k=0;
	int dist=0;
	int stringlen =  seqs[0].length();

	for (i=0; i<seqs.size(); i++){           // over row i
		for (j=0; j<=i; j++){            // over col j		
			
			for (k=0; k<stringlen; k++){ // over string length
				if (seqs[i][k] != seqs[j][k]){
					dist++;
				}
			}
			
			distmat[i][j] = dist;
			distmat[j][i] = dist;
			dist = 0;
			
		}
	}
	
}


// This function finds the minimum in a distance matrix. If there are
// ties, a pair is randomly chosen.
void findMinPair(int mypair[], vector<vector<double> >& distmat){
	
	vector<int> apair(2);
	vector<vector<int> > equal_pairs(1,apair);
	equal_pairs[0][0] = 1;
	equal_pairs[0][1] = 0;
	unsigned int i,j;
	int pairID;
	double mindist = distmat[equal_pairs[0][0]][equal_pairs[0][1]];
	
	for (i=2; i<distmat.size(); i++){ // over row
		for (j=0; j<i; j++){   // over column, but half matrix only
			
			if ( distmat[i][j] < mindist ){
				 
				// New minimum distance so replace equal_pairs with this pair
				mindist = distmat[i][j];
				equal_pairs.resize(1); 
				equal_pairs[0][0] = i;
				equal_pairs[0][1] = j;
						
			} else if ( distmat[i][j] == mindist ){ 
				
				// New pair with same distance to be added to equal_pairs
				apair[0]=i;
				apair[1]=j;
				equal_pairs.push_back(apair);
				
			}
		}
	}
	
	// Determine which pair will the minimum. If there is more than one
	// choice, randomly choose a pair.
	if (equal_pairs.size() > 1){
		//pairID = rand() % equal_pairs.size();
		random_device rd;
		pairID = rd() % equal_pairs.size();
		mypair[0] = equal_pairs[pairID][0];
		mypair[1] = equal_pairs[pairID][1];
	} else {
		mypair[0] = equal_pairs[0][0];
		mypair[1] = equal_pairs[0][1];
	}
		
	int temp;
	
	if (mypair[0]>mypair[1]){
		temp = mypair[0];
		mypair[0] = mypair[1];
		mypair[1] = temp;
	}
		
	
}


// This function is used to update the distance matrix after a 
// merge of the pair that is closest together. This assumes that 
// mypair[0] < mypair[1] which it should be by construction
void updateDistMat(vector<vector<double> >& distmat, int mypair[],
				   const vector<int>& sizes){
	
	vector<double> newdist(distmat.size()), temp_row, newdist_sub;
	vector<double>::iterator it;
	vector< vector<double> > newdistmat;
	unsigned int i,j;
	double first, second;
	
	// Create intermediate distance matrix
	for (i=0; i<newdist.size(); i++){
		first = sizes[mypair[0]]*distmat[i][mypair[0]];
		second = sizes[mypair[1]]*distmat[i][mypair[1]];
		newdist[i] = (first+second)/(sizes[mypair[0]]+sizes[mypair[1]]);
	}
	
	newdist_sub = newdist;
	it = newdist_sub.begin();
	newdist_sub.erase(it+mypair[1]);
	newdist_sub[mypair[0]] = 0;
	
	for (i=0; i<distmat.size(); i++){
		if ( (int)i!=mypair[0] && (int)i!=mypair[1] ){
			for (j=0; j<distmat.size(); j++){
				temp_row.push_back(distmat[i][j]);
			}
			it = temp_row.begin();
			temp_row.erase(it+mypair[1]);
			temp_row[mypair[0]] = newdist[i];
			newdistmat.push_back(temp_row);
			temp_row.clear();
		}
		else if ( (int)i == mypair[0] ){
			newdistmat.push_back(newdist_sub);
		} 
	}
	
	// Replace old distance matrix with the intermediate distance matrix
	distmat.pop_back(); // remove final row
	for (i=0; i<distmat.size(); i++){ // remove last column and replace contents
		distmat[i].pop_back();
		for (j=0;j<distmat[i].size();j++){
			distmat[i][j] = newdistmat[i][j];
		}
	}
	
}
