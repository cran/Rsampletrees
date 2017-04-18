/*----------------------------------------------------------------------------

  UPGMA.h - This is the header file for UPGMA.cpp

  Kelly Burkett; March 2007

  --------------------------------------------------------------------------*/

#ifndef UPGMA_H
#define UPGMA_H

#include <vector> 
#include <string>

using namespace std;

// This function sets up a distance matrix 
void computeDistMat( vector<vector<double> >& distmat, const vector<string> seqs );

// This function determines the minimum pair in a distance matrix
void findMinPair( int mypair[], vector<vector<double> >& distmat );

// This function resizes the distance matrix based on mypair
void updateDistMat( vector<vector<double> >& distmat, int mypair[], const vector<int>& sizes );

#endif
