/*----------------------------------------------------------------------------

  readFiles.h - This is the interface file for readFiles.cpp

  Kelly Burkett; March 15, 2007

  Modification Notes:

  --------------------------------------------------------------------------*/

#include <fstream>
#include <string>
#include <vector>

#include "nodefunctions.h"

#ifndef READFILES_H
#define READFILES_H



using namespace std; 

struct Options
{
	string run_name;
	unsigned long long int seed;
	int output;
	unsigned int burn_in;
	unsigned int thinning;
	unsigned long int len_chain;
	string datafile;
	int x;
	string location_file;
	int initial_tree;
	string initial_tree_file;
	int initial_haplos;
	string initial_haplo_file;
	bool random;
//	bool st_rho;
//	bool st_SA;
//	string st_file;
	char datatype;
	string haplo_freq_file;
	double initial_theta;
	double initial_rho;
	double min_theta;
	double max_theta;
	double scale_rho;
	double shape_rho;
	string weight_file;
	int haplo_list;
	string haplo_list_file;
};

struct SeqData
{
	SeqData(): seq_ID(1) {}; // don't know why I need this, something about constructor
	string seq_type;
	int num_of_type;
	vector<int> seq_ID;
};


struct genoNode
{
	string genotype;
	int id;
	int seq1, seq2;
	vector<int> hetsites;
	vector<vector<string> > hapconfigs;
};

// This function reads in the options file specified by the user
void readOptions(Options& myoptions, ifstream& myparamfile);

// This function just displays what was read in by readOptions
void printOptions(Options& myoptions, ifstream& myparamfile);

// This function reads in the sequences from a file and stores them
// in a vector of type strings
void readSeq(ifstream& myseqfile, vector<string>& seqs);

// This function reads in genotypes from a file and stores them in
// a vector of genotype data nodes. There is one node for each
// row of the data file
void readGenos( ifstream& genofile, vector<genoNode>& genodat,
			    vector<int>& hetIDs, vector<vector<int> >& hapindexmat,
			    int haploflag, const vector<string>& haplolist);

// This function reads in the distances from a file and stores them
// in a vector of distances
void readDist(ifstream& mydistfile, vector<double>& dists);

// readTree: reads in a tree in NEWICK format from a file to
// compare sampled trees to. getTimes is called by readTree and it
// extracts the times since the present from the NEWICK branch lengths
treeNode* readTree( ifstream& treefile , int id );
void getTimes ( treeNode* head_node );

// This function displays the sequences that were read in by readSeq
void printSeq(const vector<string>& seqs);

// This function returns a vector containing information about the unique
// sequences read in from the sequence file.
vector<SeqData> uniqueSeq(const vector<string> seqs);

// This function reads the haplotype frequencies from a file specified by the
// user
void readHaps( ifstream& myhapfile, vector<vector<double> >& hap_freqs,
			   unsigned int numcols);

void initialHaplos_freqs( vector<string>& seqs, vector<genoNode>& data,
					const vector<vector<double> >& hap_freqs );

void initialHaplos_enum( vector<string>& seqs, vector<genoNode>& data );

void readst( vector<double>& constants, vector<double>& ranges,
			 ifstream& stfile, Options myoption);

void readWeights( vector<double>& weights, vector<vector<int> >& updates,
		          string wfile, char type, /*bool st,*/ int& flag9 );

void findHaplos_all(vector<string>& haplos, vector<vector<string> >& allhaplos,
		string genotype, int locus, int numloci, bool firsthet );

void findHaplos_approx(bool isMissing, vector<vector<string> >& allhaps,
		genoNode sample, const vector<string>& haplolist );

#endif
