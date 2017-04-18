/*----------------------------------------------------------------------------

  p7.h -  This is the header file for the haplotype sequence update

  Kelly Burkett; September 2009

  Modification Notes:

  --------------------------------------------------------------------------*/

#ifndef P7_H
#define P7_H

#include <vector>
#include <iostream>


// The following 3 functions are called by one of the two updateP7 functions
// and they do the actual sequence change.
double updateS_haploprobs( bool sample, const vector<int>& sites,
		treeNode* chosen_ptr1, treeNode* chosen_ptr2,
		const vector<vector<double> >& haplo_probs, genoNode chosen_het );

void updateS_alleleswap( const vector<int>& siteindex, treeNode* chosen_ptr1,
		treeNode* chosen_ptr2);

int updateS_dict(vector<double>& qprobs, const vector<treeNode*>& node_ptrs,
		const vector<genoNode> samples, const double theta, vector<double>& scores,
		const vector<vector<int> >& hethaplo, vector<int>& oldindex, bool& updateflag );

int simscore(string seq, string other_seq);


// This update type changes only one locus. If updatetype = "swap" then the
// locus is changed by swapping one allele. If updatetype = "prob" then the
// new haplotypes are drawn based on the haplotype model.
int updateP7_sequence ( vector<treeNode*>& node_ptrs,
		const vector<genoNode>& samples, const vector<int>& hets,
		vector<vector<double> >& haplo_probs,
		vector<vector<double> >& allele_probs,
		double theta, vector<int>& LR_index, string updatetype, int nloci,
		bool& updateflag, vector<string>& newseqs);


double probP7_sequence( genoNode chosen_het, treeNode* node_ptr_1, treeNode* node_ptr_2,
		treeNode* copy_node_1, treeNode* copy_node_2, double theta,
		vector<int>& first_markers, vector<vector<double> >& allele_probs,
		vector<vector<double> >& hap_probs );


// This update is like updateP7_sequence except that the nodes chosen to be
// updated are sampled with unequal weights
int updateP7_weighted ( vector<treeNode*>& node_ptrs, const vector<genoNode>& samples,
						const vector<int>& hets,
						vector<vector<double> >& haplo_probs,
						vector<vector<double> >& allele_probs,
						double theta, vector<int>& LR_index,
						vector<double>& het_weights, string updatetype, int nloci );


// This performs the sequence update plus topology change. The topology change
// is done with updates to the r/s values up the tree.
int updateP7_topo ( vector<treeNode*>& node_ptrs,
					 vector<genoNode>& samples,
			   		 vector<int>& hets,
	  		   		 vector<vector<double> >& haplo_probs,
   			   		 vector<vector<double> >& allele_probs,
   			   		 double theta, double rho,
					 vector<int>& first_markers,
					 vector<double>& locations, int ref_marker,
					 vector<double>& scores,
					 const vector<vector<int> >& hethaplo,
					 vector<int>& oldindex,
					 bool& updateflag );


void selectOtherNodeP7( vector<double>& probs, vector<treeNode*>& values,
						vector<treeNode*>& node_ptrs, treeNode* chosen_node,
					    vector<int>& first_markers, double theta, double rho,
						   const vector<double>& locations, treeNode* other_node=NULL );

double updateUpTree( bool sample, treeNode* chosen_1, treeNode* chosen_2,
					 treeNode* oldsib_1, treeNode* oldsib_2,
	  				 treeNode* other_1, treeNode* other2,
			  	     const double rho, const double theta,
   					 vector<int>& first_markers,
   					 vector<double>& locations, int ref_marker,
   					 vector<vector<double> >& allele_probs,
   					 vector<vector<double> >& hap_probs);

void topologyChangeP7( treeNode* node_ptr, treeNode* other_ptr );

double probP7_major( treeNode* node_ptr1, treeNode* other_ptr1, treeNode* sib_ptr1,
					 treeNode* node_ptr2, treeNode* other_ptr2, treeNode* sib_ptr2,
			   	     vector<treeNode*> node_ptrs, double theta, double rho,
   					 vector<int>& first_markers, int ref_marker,
   					 vector<double>& locations, vector<vector<double> >& allele_probs,
   					 vector<vector<double> >& hap_probs,
   					 int min_index, int max_index );

double updateS_P7(bool sample, treeNode* node_ptr,
				  vector<int>& first_markers, double theta,
	  			  vector<vector<double> >& allele_probs,
				  vector<vector<double> >& hap_probs);

void getSijPMF_P7(string side, treeNode* node_ptr, vector<double>& probs,
			   vector<int>& first_markers, int locus, double theta,
	  		   vector<vector<double> >& allele_probs,
   			   vector<vector<double> >& hap_probs);

double probSij_P7(string side, treeNode* node_ptr, int locus, double theta,
			   vector<int>& first_markers, vector<vector<double> >& allele_probs,
	 		   vector<vector<double> >& hap_probs);


void swapNodes(genoNode chosen_node, vector<treeNode*>& node_ptrs);

void updateHaploList( const vector<string>& newseqs, vector<string>& haplo_list_vector,
		vector<genoNode>& samples, vector< vector<int> >& hapindexmat,
		vector<double>& P9scores );

#endif
