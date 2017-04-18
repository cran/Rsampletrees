/*----------------------------------------------------------------------------

  p7.cpp - This file contains functions needed to complete the haplotype
  	       sequence update

  Kelly Burkett; September 2009

  Modification Notes:

  --------------------------------------------------------------------------*/

#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <bitset>
#include <gsl/gsl_randist.h>
#include <Rcpp.h>

#include "nodefunctions.h"
#include "rng.h"
#include "readFiles.h"
#include "models.h"
#include "p7.h"
#include "misc.h"
#include "proposal.h"
#include "treeBuild.h"


// *************************************************************************//
// SEQUENCE UPDATE FUNCTIONS USED BY ALL HAPLOTYPE UPDATE TYPES
// *************************************************************************//

// This function updates the sequence using an allele swap at the loci in
// the vector sites. Since the swap is symmetric, the function doesn't need
// to return the probability
void updateS_alleleswap( const vector<int>& sites, treeNode* chosen_ptr1,
		treeNode* chosen_ptr2)
{

	int i=0, numsites=sites.size();
	vector<double> probs(2,0.5);
	int change=1;

	for (i=0; i<numsites; i++){

		//change=samplePMF( probs );

		if (change==1){
			if ( chosen_ptr1->seq[sites[i]] == '0' ){
				if (chosen_ptr2->seq[sites[i]] != '1' ){
					/*cerr*/ Rcpp::Rcout << "ERROR - A homozygote was sampled in the haplotype sampler" << endl;
					throw Rcpp::exception("ERROR - A homozygote was sampled in the haplotype sampler"); //exit(1);
				} else {
					chosen_ptr1->seq[sites[i]] = '1';
					chosen_ptr2->seq[sites[i]] = '0';
				}
			} else {
				if (chosen_ptr2->seq[sites[i]] != '0' ){
					/*cerr*/ Rcpp::Rcout << "ERROR - A homozygote was sampled in the haplotype sampler" << endl;
					throw Rcpp::exception("ERROR - A homozygote was sampled in the haplotype sampler"); //exit(1);
				} else {
					chosen_ptr1->seq[sites[i]] = '0';
					chosen_ptr2->seq[sites[i]] = '1';
				}
			}
		}

	}
}


// This function updates the sequence by choosing the haplotype phase from the
// haplotype probability model for each of the loci in sites.
// Since this update is not symmetric, it returns the probability of the
// change and can be used to find the probability of the current values
// if sample=0
double updateS_haploprobs( bool sample, const vector<int>& sites,
		treeNode* chosen_ptr1, treeNode* chosen_ptr2,
		const vector<vector<double> >& haplo_probs, genoNode chosen_het )
{

	int i=0, j=0, numsites=sites.size();
	double prob=1, prob00=0, unif00=0;
	int allele=0;

	// Loop to change the sequence or determine the probabilities of the
	// current sequence
	for (i=0; i<numsites; i++){

		j = sites[i];

		if ( (j==0) || (chosen_het.genotype[j-1]=='0') ||
					(chosen_het.genotype[j-1]=='2') ){

			// Only one two locus haplotype configuration in these cases. Have
			// equal probabilities.
			if ( sample==0 ){
				prob = prob * 0.5;
			} else {
				allele = generateDiscreteUnif(0,1);
				chosen_ptr1->seq[j]=itoc(allele);
				chosen_ptr2->seq[j]=itoc(1-allele);
				prob = prob * 0.5;
			}

		} else {

			// Two consecutive heterozygous SNPs. Sample new configuration
			// based on the haplotype probabilities and the weights of the
			// two configurations


			// Determine the probabilities of the two configurations
			// prob00 is Pr(00/11)/[Pr(00/11)+Pr(01/10)].
			prob00 = (haplo_probs[j-1][0]*haplo_probs[j-1][3])/
					 (haplo_probs[j-1][0]*haplo_probs[j-1][3]+
					  haplo_probs[j-1][1]*haplo_probs[j-1][2]);


			if ( sample == 0 ){ // Find probability of current sequence

				if ( chosen_ptr1->seq[j-1]=='1' ){
					if (chosen_ptr1->seq[j] == '1'){// seq1 is 11 so seq2 is 00
						prob = prob * prob00;
					} else {             // seq1 is 10 so seq2 is 01
						prob = prob * (1-prob00);
					}
				} else {
					if (chosen_ptr1->seq[j] == '0'){// seq1 is 00 so seq1 is 11
						prob = prob * prob00;
					} else {             // seq1 is 01 so seq1 is 10
						prob = prob * (1-prob00);
					}
				}
			} else { // Propose new sequence data

				unif00 = generateStdUnif();

				if (unif00<prob00){ // Configuration selected is 00/11
					if ( chosen_ptr1->seq[j-1]=='1' ){// seq1 is 11
						chosen_ptr1->seq[j]='1';
						chosen_ptr2->seq[j]='0';
					} else {                          // seq1 is 00
						chosen_ptr1->seq[j]='0';
						chosen_ptr2->seq[j]='1';
					}
					prob = prob * prob00;
				} else { // Configuration selected is 01/10
					if ( chosen_ptr1->seq[j-1]=='1' ){// seq1 is 10
						chosen_ptr1->seq[j]='0';
						chosen_ptr2->seq[j]='1';
					} else {                          // seq1 is 01
						chosen_ptr1->seq[j]='1';
						chosen_ptr2->seq[j]='0';
					}
					prob = prob * (1-prob00);
				}
			}
		}


	}

	return(prob);

}


int updateS_dict(vector<double>& qprobs, const vector<treeNode*>& node_ptrs,
		const vector<genoNode> samples, const double theta, vector<double>& scores,
		const vector<vector<int> >& hethaplo, vector<int>& oldindex, bool& updateflag )
{

	int numhaplos=0;
	int numsamples=samples.size();
	int numtips=(node_ptrs.size()+1)/2;
	int seqlength=node_ptrs[0]->seq.size();
	int i=0, j=0, k=0, index=0, max1=0, max2=0,stop=0;
	string seq1="", seq2="", other_seq1="", other_seq2="";
	genoNode chosen_het;
	int chosen_haps, chosenid;
	int counter = 0;
	vector<double> probs(scores.size(),0);
	vector<double> scaledprobs(scores.size(),0);
	double total=0;


	if (updateflag == 1){ // Recompute the scores and the total

		for (i=0; i<numsamples; i++){

			if ( samples[i].hetsites.size()>1 ){

				seq1 = node_ptrs[samples[i].seq1]->seq;
				seq2 = node_ptrs[samples[i].seq2]->seq;
				numhaplos = samples[i].hapconfigs.size();

				for (k=0; k<numhaplos; k++){

					other_seq1 = samples[i].hapconfigs[k][0];
					other_seq2 = samples[i].hapconfigs[k][1];

					if ( (seq1==other_seq1)||(seq1==other_seq2) ){
						// The potential configuration is already the current configuration
						// Make this have probability 0 of being chosen.
						scores[counter]=0;
						oldindex[i]=counter;
					} else {
						j=0;
						while( (j<numtips)&&(stop==0) ){
							// Compute the similarity score
							max1 = max(max1,simscore(other_seq1,node_ptrs[j]->seq));
							max2 = max(max2,simscore(other_seq2,node_ptrs[j]->seq));
							if ( (max1==seqlength)&&(max2==seqlength) ){
								stop=1;
							}
							j++;
						}

						scores[counter]=max1+max2;
						max1=0;
						max2=0;
						stop=0;

					}
					counter++;

				}

			} else {
				oldindex[i]=-1;
				scores[counter]=0;
				counter++;
			}


		}

	}


	// Compute the similarity scores. They depend on theta
	// Rescale the similarity score
	for (i=0; i<(int)scores.size(); i++){
		if (scores[i]!=0){
			probs[i]=gsl_ran_poisson_pdf(2*seqlength-scores[i],theta);
			total = total+probs[i];
		} else {
			probs[i] = 0;
		}
	}

	for (i=0; i<(int)scores.size(); i++){
		scaledprobs[i]=probs[i]/total;
	}


	// Sample a new configuration and determine the probability
	index = samplePMF(scaledprobs);
	qprobs[1] = scaledprobs[index];
	chosenid = hethaplo[index][0];
	chosen_het = samples[chosenid];

	if (chosen_het.hetsites.size()==0){
		/*cerr*/ Rcpp::Rcout << "A homozygote was chosen" << endl;
		throw Rcpp::exception("A homozygote was chosen"); //exit(1);
	}
	chosen_haps = hethaplo[index][1];
	node_ptrs[chosen_het.seq1]->seq = chosen_het.hapconfigs[chosen_haps][0];
	node_ptrs[chosen_het.seq2]->seq = chosen_het.hapconfigs[chosen_haps][1];


	// Find the probability of the old sequence. For this we need to
	// compute the probability for the old configuration since it
	// was previously 0, and remove the probability of the new configuration
	double oldscore=0;

	j=0;
	chosen_haps=hethaplo[oldindex[chosenid]][1];
	other_seq1 = chosen_het.hapconfigs[chosen_haps][0];
	other_seq2 = chosen_het.hapconfigs[chosen_haps][1];
	while( (j<numtips)&&(stop==0) ){
		// Compute the similarity score
		max1 = max(max1,simscore(other_seq1,node_ptrs[j]->seq));
		max2 = max(max2,simscore(other_seq2,node_ptrs[j]->seq));
		if ( (max1==seqlength)&&(max2==seqlength) ){
			stop=1;
		}
		j++;

	}


	// Remove term corresponding to the proposed sequence from the total
	// Multiply by the old total to get the unscaled value.
	// In A~->A this term would be 0 but I don't want to recompute everything.
	oldscore = gsl_ran_poisson_pdf(2*seqlength-max1-max2,theta);


	// Now get score for A~->A by getting the new standardizing constant
	// and dividing the score by that value.
	total = total - probs[index] + oldscore;
	qprobs[0] = oldscore/total;

	return(chosen_het.id);

}


int simscore(string seq, string other_seq)
{
	int length=seq.size();
	int score=0;

	for (int i=0; i<length; i++){
		if (seq[i]==other_seq[i]){
			score+=1;
		}
	}
	return(score);
}


// *************************************************************************//
// HAPLOTYPE UPDATE / NO TOPO CHANGE FUNCTIONS
// *************************************************************************//

// This is an update function that does not make a topology change. It updates
// only the sequence of a randomly chosen phase-unknown individual. The sequence
// update is made to nloci of the heterozygous sites and the sequence change
// is chosen based on the value of updatetype.
// Return the new sequence just in case the dictionary rephaser is being used
// as well.
int updateP7_sequence ( vector<treeNode*>& node_ptrs,
		const vector<genoNode>& samples, const vector<int>& hets,
		vector<vector<double> >& haplo_probs,
		vector<vector<double> >& allele_probs, double theta,
		vector<int>& LR_index, string updatetype, int nloci, bool& updateflag,
		vector<string>& newseqs)
{
	int accept = 1;
	double prob=log(double(1)), MHratio=1, qA=1, qAnew=1;
	double myunif = generateStdUnif();


	// Choose an element from the vector of unphased individual labels.
	// Then extract information for that individual from the vector of
	// genotype info.
	int hetindex = generateDiscreteUnif( 0,hets.size()-1 ) ;
	genoNode chosen_het = samples [ hets[hetindex]-1 ];


	// Determine node labels and the sequences for the two nodes corresponding
	// to the chosen individuals two sequences
	int seq1index = chosen_het.seq1;
	int seq2index = chosen_het.seq2;
	treeNode copy_node_1 = *(node_ptrs[seq1index]);
	treeNode copy_node_2 = *(node_ptrs[seq2index]);


	// Sample nloci heterozygous sites in this individual. Store it
	// as a vector because there different ways to complete this update
	vector<int> copy_hetsites = chosen_het.hetsites;
	vector<int> sites(0);
	int siteindex=-1;

	if ( nloci>=(int)chosen_het.hetsites.size() ){
		sites=chosen_het.hetsites;
	} else {
		for (int i=0; i<nloci; i++){
			siteindex = generateDiscreteUnif( 0, copy_hetsites.size()-1 );
			sites.push_back(copy_hetsites[siteindex]);
			copy_hetsites.erase(copy_hetsites.begin()+siteindex);
		}
		sort (sites.begin(), sites.end() );
	}


	// Make the change to the sequence
	if (updatetype == "swap"){
		updateS_alleleswap(sites, node_ptrs[seq1index],
				node_ptrs[seq2index]);
	} if (updatetype == "prob"){
		qA = updateS_haploprobs( 0, sites, node_ptrs[seq1index],
				node_ptrs[seq2index], haplo_probs, chosen_het );
		qAnew = updateS_haploprobs( 1, sites, node_ptrs[seq1index],
				node_ptrs[seq2index], haplo_probs,chosen_het );
	}


	// Compute log(prob(A~|G)/Pr(A|G)) for this update
	prob = probP7_sequence( chosen_het, node_ptrs[seq1index], node_ptrs[seq2index],
			&copy_node_1, &copy_node_2, theta, LR_index,
			allele_probs, haplo_probs );


	// Compute the acceptance probability and determine whether the change is
	// accepted or not. If not, replace the new sequences with the old.
	MHratio = min( exp(prob)*qA/qAnew, double( 1 ) );

	if ( myunif > MHratio ){
		accept=0;
		node_ptrs[seq1index]->seq = copy_node_1.seq;
		node_ptrs[seq2index]->seq = copy_node_2.seq;
	} else {
		// Test to see if the node labels should be changed. They are changed
		// if the update to the sequence has made it so that the nodes need
		// to be relabelled. If they need to be relabelled then we also need
		// to change node_ptrs and the sequences
		if (node_ptrs[seq1index]->seq[chosen_het.hetsites[0]]=='0'){
			swapNodes(chosen_het, node_ptrs);
		}
		updateflag=1;

		// If an update is accepted, return the new sequences just in case
		// the dictionary rephaser is being used.
		newseqs[0]=node_ptrs[seq1index]->seq;
		newseqs[1]=node_ptrs[seq2index]->seq;
	}

	return( accept );

}


// This computes the term log(Pr(A~|G)/Pr(A|G)) = log Pr(A~|G) - log Pr(A|G)
// All terms are initially stored in a vector for the numerator and
// denominator, then each vector is sorted according to size and the smallest
// elements are first subtracted from each other to help ensure that there
// aren't any numerical errors
double probP7_sequence( genoNode chosen_het, treeNode* node_ptr_1, treeNode* node_ptr_2,
		treeNode* copy_node_1, treeNode* copy_node_2, double theta,
		vector<int>& first_markers, vector<vector<double> >& allele_probs,
		vector<vector<double> >& hap_probs )
{

	double logprobRatio = log(double(1));
	vector<double> probsA(0);
	vector<double> probsAnew(0);
	int i=0, nsites=chosen_het.hetsites.size();
	string side="left";
 	int nloci = node_ptr_1->seq.size();

	for (int j=0; j<nsites; j++){

		i=chosen_het.hetsites[j];

		if (i >= first_markers[1]){
			side="right";
		}

		// Compute the probability for the locus that has changed
		probsA.push_back(log(probSij( side, copy_node_1, i, theta, first_markers, allele_probs, hap_probs)));
		probsA.push_back(log(probSij( side, copy_node_2, i, theta, first_markers, allele_probs, hap_probs)));

		probsAnew.push_back(log(probSij( side, node_ptr_1, i, theta, first_markers, allele_probs, hap_probs)));
		probsAnew.push_back(log(probSij( side, node_ptr_2, i, theta, first_markers, allele_probs, hap_probs)));

		if (side=="left"){

			// Compute the probability for the "next" locus (since
			// side="left" the next is actually to the left) since
			// if its prob comes from the haplo model its value
			// will change too
			if ( ((i-1)!=-1)&&((i-1)!=chosen_het.hetsites[j-1]) ){

				probsA.push_back(log(probSij( side, copy_node_1, i-1, theta, first_markers, allele_probs, hap_probs)));
				probsA.push_back(log(probSij( side, copy_node_2, i-1, theta, first_markers, allele_probs, hap_probs)));

				probsAnew.push_back(log(probSij( side, node_ptr_1, i-1, theta, first_markers, allele_probs, hap_probs)));
				probsAnew.push_back(log(probSij( side, node_ptr_2, i-1, theta, first_markers, allele_probs, hap_probs)));
			}
		} else {

			// Compute the probability for the "next" locus (since
			// side="right" the next is to the right) since
			// if its prob comes from the haplo model its value
			// will change too
			if ( ((i+1)!=nloci)&&((i+1)!=chosen_het.hetsites[j+1]) ){

				probsA.push_back(log(probSij( side, copy_node_1, i+1, theta, first_markers, allele_probs, hap_probs)));
				probsA.push_back(log(probSij( side, copy_node_2, i+1, theta, first_markers, allele_probs, hap_probs)));

				probsAnew.push_back(log(probSij( side, node_ptr_1, i+1, theta, first_markers, allele_probs, hap_probs)));
				probsAnew.push_back(log(probSij( side, node_ptr_2, i+1, theta, first_markers, allele_probs, hap_probs)));
			}
		}


	}

	// Sort both vectors from smallest to largest
	sort(probsA.begin(), probsA.end());
	sort(probsAnew.begin(), probsAnew.end());

	// Compute the log ratio by adding smallest differences first
	int nterms = probsA.size();

	for (i=0; i<nterms; i++){
		logprobRatio = logprobRatio+probsAnew[i]-probsA[i];
	}

	return( logprobRatio );

}






int updateP7_weighted ( vector<treeNode*>& node_ptrs, const vector<genoNode>& samples,
						const vector<int>& hets,
						vector<vector<double> >& haplo_probs,
						vector<vector<double> >& allele_probs,
						double theta, vector<int>& LR_index,
						vector<double>& het_weights, string updatetype, int nloci )
{
	int accept = 1;
	double MHratio=1, probRatio=log(double(1)), qA=1, qAnew=1;
	double myunif = generateStdUnif();
	vector<double> probsA;
	vector<double> probsAnew;


	int hetindex = samplePMF( het_weights ) ;
	genoNode chosen_het = samples [ hetindex ];


	int seq1index = chosen_het.seq1;
	int seq2index = chosen_het.seq2;
	string seq1 = node_ptrs[seq1index]->seq;
 	string seq2 = node_ptrs[seq2index]->seq;
 	int nsites = seq1.size();

 	string side="left";
 	treeNode copy_node1= *(node_ptrs[seq1index]);
 	treeNode copy_node2= *(node_ptrs[seq2index]);


 	// Sample nloci heterozygous sites in this individual. Store it
 	// as a vector because there different ways to complete this update
 	vector<int> copy_hetsites = chosen_het.hetsites;
 	vector<int> sites(0);
 	int siteindex=-1;

 	if ( nloci>=(int)chosen_het.hetsites.size() ){
 		sites=chosen_het.hetsites;
	} else {
		for (int i=0; i<nloci; i++){
			siteindex = generateDiscreteUnif( 0, copy_hetsites.size()-1 );
			sites.push_back(copy_hetsites[siteindex]);
			copy_hetsites.erase(copy_hetsites.begin()+siteindex);
		}
		sort (sites.begin(), sites.end() );
	}


 	// Make the change to the sequence
 	if (updatetype == "swap"){
 		updateS_alleleswap(sites, node_ptrs[seq1index],
 				node_ptrs[seq2index]);
 	} if (updatetype == "prob"){
 		qA = updateS_haploprobs( 0, sites, node_ptrs[seq1index],
 				node_ptrs[seq2index], haplo_probs, chosen_het );
 		qAnew = updateS_haploprobs( 1, sites, node_ptrs[seq1index],
 				node_ptrs[seq2index], haplo_probs,chosen_het );
 	}


	// Compute the probability of the change in a safe manner
	int i=0;

	for (int j=0; j<(int)chosen_het.hetsites.size(); j++){

		i=chosen_het.hetsites[j];

		if (i >= LR_index[1]){
			side="right";
		}

		probsA.push_back(log(probSij( side, &copy_node1, i, theta, LR_index, allele_probs, haplo_probs)));
		probsA.push_back(log(probSij( side, &copy_node2, i, theta, LR_index, allele_probs, haplo_probs)));

		probsAnew.push_back(log(probSij( side, node_ptrs[seq1index], i, theta, LR_index, allele_probs, haplo_probs)));
		probsAnew.push_back(log(probSij( side, node_ptrs[seq2index], i, theta, LR_index, allele_probs, haplo_probs)));

		if (side=="left"){
			if ( ((i-1)!=-1)&&((i-1)!=chosen_het.hetsites[j-1]) ){

				probsA.push_back(log(probSij( side, &copy_node1, i-1, theta, LR_index, allele_probs, haplo_probs)));
				probsA.push_back(log(probSij( side, &copy_node2, i-1, theta, LR_index, allele_probs, haplo_probs)));

				probsAnew.push_back(log(probSij( side, node_ptrs[seq1index], i-1, theta, LR_index, allele_probs, haplo_probs)));
				probsAnew.push_back(log(probSij( side, node_ptrs[seq2index], i-1, theta, LR_index, allele_probs, haplo_probs)));
			}
		} else {
			if ( ((i+1)!=nsites)&&((i+1)!=chosen_het.hetsites[j+1]) ){

				probsA.push_back(log(probSij( side, &copy_node1, i+1, theta, LR_index, allele_probs, haplo_probs)));
				probsA.push_back(log(probSij( side, &copy_node2, i+1, theta, LR_index, allele_probs, haplo_probs)));

				probsAnew.push_back(log(probSij( side, node_ptrs[seq1index], i+1, theta, LR_index, allele_probs, haplo_probs)));
				probsAnew.push_back(log(probSij( side, node_ptrs[seq2index], i+1, theta, LR_index, allele_probs, haplo_probs)));
			}
		}


	}

	sort(probsA.begin(), probsA.end());
	sort(probsAnew.begin(), probsAnew.end());

	for (int j=0; j<(int)probsA.size(); j++){
		probRatio = probRatio+probsAnew[j]-probsA[j];
	}

	MHratio = min( exp(probRatio)*qA/qAnew, double(1) );


	if ( myunif > MHratio ){
		accept=0;
		node_ptrs[seq1index]->seq = seq1;
		node_ptrs[seq2index]->seq = seq2;
	} else {
		accept=0;
		if ( (node_ptrs[seq1index]->seq != seq1)&&(node_ptrs[seq1index]->seq != seq2) ){
			accept=1;
		}
	}


	return( accept );
}





// *************************************************************************//
// HAPLOTYPE UPDATE / WITH TOPO CHANGE FUNCTIONS
// *************************************************************************//
int updateP7_topo ( vector<treeNode*>& node_ptrs,
					 vector<genoNode>& samples,
			   		 vector<int>& hets,
	  		   		 vector<vector<double> >& haplo_probs,
   			   		 vector<vector<double> >& allele_probs,
   			   		 double theta, double rho,
					 vector<int>& first_markers,
					 vector<double>& locations, int ref_marker,
					 vector<double>& scores,
					 const vector<vector<int> >& hethaplo, vector<int>& oldindex,
					 bool& updateflag )
{
	int accept = 1, index=0;
	double probA=1, probAnew=1, qA=1, qAnew=1, MHratio=1;
	double myunif = generateStdUnif(), temp=0;
	string seq1="", seq2="";
	int n = (node_ptrs.size()+1)/2;


	// Make a copy of the tree to revert back to if the update is not accepted
	// and to use for computing the probability of the tree before the change
	treeNode* copy_head_node = copyTree(node_ptrs[node_ptrs.size()-1]);
	vector<treeNode*> copy_node_ptrs(node_ptrs.size());
	getNodePtrs(copy_node_ptrs, copy_head_node);

	treeNode* reverse_copy_head_node = copyTree(node_ptrs[node_ptrs.size()-1]);
	vector<treeNode*> reverse_copy_node_ptrs(node_ptrs.size());
	getNodePtrs(reverse_copy_node_ptrs, reverse_copy_head_node);


	// Propose new sequence for the two nodes. Find the probability
	// associated with the original sequence as well.
	vector<double> qprobs(2,0);
	int chosen_index=-1;

	chosen_index=updateS_dict(  qprobs, node_ptrs, samples, theta,
			scores, hethaplo, oldindex, updateflag )-1;

	if (updateflag==1){
		// We just updated the scores, so they won't need to be updated on the
		// next round
		updateflag=0;
	}
	qAnew = qAnew * qprobs[1];
	qA = qA * qprobs[0];


	// Setup for the individual chosen
	treeNode *chosen_ptr1 = node_ptrs[samples[chosen_index].seq1];
	treeNode *chosen_ptr2 = node_ptrs[samples[chosen_index].seq2];


	// Get pointers to the chosen nodes in the copy of the tree. Since
	// chosen_het.seq gives the index in the vector node_ptrs, it also
	// gives the index in the vector copy_node_ptrs
	treeNode *copy_chosen_ptr1 = copy_node_ptrs[samples[chosen_index].seq1];
	treeNode *copy_chosen_ptr2 = copy_node_ptrs[samples[chosen_index].seq2];

	treeNode *reverse_copy_chosen_ptr1 = reverse_copy_node_ptrs[samples[chosen_index].seq1];
	treeNode *reverse_copy_chosen_ptr2 = reverse_copy_node_ptrs[samples[chosen_index].seq2];


	// Determine the sibs of the chosen nodes in the initial topology
	// and in its copy. The copy will remain unchanged.
	treeNode *sib_ptr1 = NULL, *copy_sib_ptr1=NULL, *reverse_copy_sib_ptr1=NULL;
	treeNode *sib_ptr2 = NULL, *copy_sib_ptr2=NULL, *reverse_copy_sib_ptr2=NULL;

	if (chosen_ptr1->parent->child1->id == chosen_ptr1->id){
		sib_ptr1 = chosen_ptr1->parent->child2;
		copy_sib_ptr1 = copy_chosen_ptr1->parent->child2;
		reverse_copy_sib_ptr1 = reverse_copy_chosen_ptr1->parent->child2;
	} else {
		sib_ptr1 = chosen_ptr1->parent->child1;
		copy_sib_ptr1 = copy_chosen_ptr1->parent->child1;
		reverse_copy_sib_ptr1 = reverse_copy_chosen_ptr1->parent->child1;
	}

	if (chosen_ptr2->parent->child1->id == chosen_ptr2->id){
		sib_ptr2 = chosen_ptr2->parent->child2;
		copy_sib_ptr2 = copy_chosen_ptr2->parent->child2;
		reverse_copy_sib_ptr2 = reverse_copy_chosen_ptr2->parent->child2;
	} else {
		sib_ptr2 = chosen_ptr2->parent->child1;
		copy_sib_ptr2 = copy_chosen_ptr2->parent->child1;
		reverse_copy_sib_ptr2 = reverse_copy_chosen_ptr2->parent->child1;
	}


	// Determine if the two nodes to be moved are related in some way. If they
	// are it might be a special case for the reverse topology rearrangement
	bool /*sibs=0,*/ niece1=0/*, aunt1=0*/;

	if ( sib_ptr1 == chosen_ptr2 ){
		// The two nodes to be moved are siblings.

//		sibs = 1;

		if (chosen_ptr1->parent->parent->child1 == chosen_ptr1->parent){

			// In order to compute prAgivenG and updateUpTree for the
			// reverse rearrangement, we need to modify the sib nodes so that
			// the correct coalescent events occur. In going from A~ -> A
			// chosen_ptr1 must connect with its old aunt, so we must replace
			// the sib_ptr1 copies with the aunt node
			copy_sib_ptr1 = copy_chosen_ptr1->parent->parent->child2;
			reverse_copy_sib_ptr1 = reverse_copy_chosen_ptr1->parent->parent->child2;

			// After the first topology change, the sib of chosen_ptr2 will
			// in fact be its old aunt. In order to ensure that the updateUpTree is
			// done to the correct branches, we need to ensure that the old aunt
			// is in the path by updating the sibptr node
			sib_ptr2 = chosen_ptr1->parent->parent->child2;


		} else {
			copy_sib_ptr1 = copy_chosen_ptr1->parent->parent->child1;
			reverse_copy_sib_ptr1 = reverse_copy_chosen_ptr1->parent->parent->child1;
			sib_ptr2 = chosen_ptr1->parent->parent->child1;
		}
	}

	if ( sib_ptr1 == chosen_ptr2->parent ){
		// The first one to move is the aunt of the second. After the topology
		// rearrangement sib_ptr1 will have moved along with chosen_ptr2 as it
		// is the parent of chosen_ptr2. Therefore, for the reverse
		// rearrangement, the copy needs to be the NIECE of chosen_ptr1, which
		// is also the sib of chosen_ptr2
		// NOTE: This was pointing to sib_ptr2... rather than copy_sib_ptr2.
		//       I think it should be pointing to the copy not the original
//		aunt1 = 1;
		copy_sib_ptr1 = copy_sib_ptr2;
		reverse_copy_sib_ptr1 = reverse_copy_sib_ptr2;
	}

	if ( sib_ptr2 == chosen_ptr1->parent ){
		// The first one to move is the niece of the second. This is only a
		// concern if they end up as sibs in the new topology so for now just
		// flag it
		niece1 = 1;
	}




	// Given the two changed sequences at chosen_ptr1 and chosen_ptr2, begin
	// the update to topology and associated variables. A step is completed
	// first for chosen_ptr1, then for chosen_ptr2, with probabilities of
	// reverse rearrangement computed where required by the intermediate changes.
	vector<double> probs(0);
	vector<treeNode*> nodes(0);
	treeNode *other_ptr1 = NULL;
	treeNode *other_ptr2 = NULL;
	treeNode *copy_other_ptr1 = NULL;
	treeNode *copy_other_ptr2 = NULL;
	treeNode *reverse_copy_other_ptr1 = NULL;
	treeNode *reverse_copy_other_ptr2 = NULL;
	treeNode *max_ptr = node_ptrs[node_ptrs.size()-1] ;
	vector<treeNode*>::iterator it;
	bool move_head = 0;

	// 1. (a) Select a new node to coalesce with chosen_ptr1.
	selectOtherNodeP7( probs, nodes, node_ptrs, chosen_ptr1,
					   first_markers, theta, rho, locations, chosen_ptr2 );
	if ( probs.size()==0 ){
		// This update can't be done so return
		return( -1 );
	}
	index = samplePMF( probs );
	other_ptr1 = nodes[index];
	copy_other_ptr1 = copy_node_ptrs[other_ptr1->id-1];
	reverse_copy_other_ptr1 = reverse_copy_node_ptrs[other_ptr1->id-1];
	qAnew = qAnew * probs[index];

	if ( (chosen_ptr1->parent->parent!=NULL)&&(other_ptr1->parent!=NULL) ){
		if (chosen_ptr1->parent->parent->t > other_ptr1->parent->t){
			max_ptr = chosen_ptr1->parent->parent;
		} else { max_ptr= other_ptr1->parent; }
	} else {
		move_head = 1;
	}

	// 1. (b) Make the topology change and update the time
	topologyChangeP7(chosen_ptr1, other_ptr1);
	topologyChangeP7(reverse_copy_chosen_ptr1, reverse_copy_other_ptr1);
	temp = updateT( 1, chosen_ptr1->parent );
	qAnew = qAnew * temp;
	changeT(reverse_copy_chosen_ptr1->parent, chosen_ptr1->parent->t);


	// 2. (a) Select a new node to coalesce with chosen_ptr2
	selectOtherNodeP7( probs, nodes, node_ptrs, chosen_ptr2,
					   first_markers, theta, rho, locations );
	if ( probs.size()==0 ){
		// This update can't be done so return
		return( -1 );
	}
	index = samplePMF( probs );
	other_ptr2 = nodes[index];
	copy_other_ptr2 = copy_node_ptrs[other_ptr2->id-1];
	reverse_copy_other_ptr2 = reverse_copy_node_ptrs[other_ptr2->id-1];
	qAnew = qAnew * probs[index];

	if ( (chosen_ptr2->parent->parent!=NULL)&&(other_ptr2->parent!=NULL) ){
		if (chosen_ptr2->parent->parent->t > other_ptr2->parent->t){
			if ( chosen_ptr2->parent->parent->t > max_ptr->t ){
				max_ptr = chosen_ptr2->parent->parent;
			}
		} else {
			if (other_ptr2->parent->t > max_ptr->t){
				max_ptr= other_ptr2->parent;
			}
		}
	} else {
		move_head = 1;
	}


	// 2. (b) Make the topology change and update the time
	topologyChangeP7(chosen_ptr2, other_ptr2);
	topologyChangeP7(reverse_copy_chosen_ptr2, reverse_copy_other_ptr2);
	temp = updateT( 1, chosen_ptr2->parent );
	qAnew = qAnew * temp;
	changeT(reverse_copy_chosen_ptr2->parent, chosen_ptr2->parent->t);



	// If chosen_ptr1 is the niece of chosen_ptr2 and they end up as sibs in
	// A~ then the parents of these two nodes will swap. In the reverse
	// rearrangement the second then needs to coalesce with the new parent
	// of chosen_ptr1
	if ( (other_ptr2 == chosen_ptr1)&&(niece1==1) ){
		reverse_copy_sib_ptr2 = reverse_copy_chosen_ptr1->parent;
	}


	// 3. Update the r and s variables for the 4 lineages up the tree that may
	//    have changed
	temp = updateUpTree( 1, chosen_ptr1, chosen_ptr2,
			  other_ptr1, other_ptr2,
			  sib_ptr1, sib_ptr2,
			  rho, theta, first_markers,
			  locations, ref_marker, allele_probs,
			  haplo_probs );
	qAnew = qAnew * temp;


	// 6. Now mirror steps 1. and 2. with the copy of the tree to find Q(A|A~)

	// Given A~, determine the probability of selecting the old sib (sib_ptr1)
	// as the sibling of the first chosen node (chosen_ptr1). This must be
	// calculated after the two topology changes have been made because we
	// want A~->A

	// NOTE: I think these lines are not necessary because we don't change the sequence in the reverse copy
	if (reverse_copy_chosen_ptr1->seq != copy_chosen_ptr1->seq){
		/*cout*/ Rcpp::Rcout << "\n\n\nNOT THE SAME!!!!!\n\n\n" << endl;
	}
	reverse_copy_chosen_ptr1->seq = copy_chosen_ptr1->seq;
	reverse_copy_chosen_ptr2->seq = copy_chosen_ptr2->seq;

	selectOtherNodeP7( probs, nodes, reverse_copy_node_ptrs, reverse_copy_chosen_ptr1,
 					   first_markers, theta, rho, locations, reverse_copy_chosen_ptr2 );
	it = find( nodes.begin(), nodes.end(), reverse_copy_sib_ptr1);
	if ( it == nodes.end() ){
		/*cerr*/ Rcpp::Rcout << "ERROR: Reverse change is not possible in Major Haplotype update" << endl;
		throw Rcpp::exception("ERROR: Reverse change is not possible in Major Haplotype update"); //exit(1);
	}
	qA = qA * probs.at(int( it - nodes.begin() ));

	topologyChangeP7(reverse_copy_chosen_ptr1, reverse_copy_sib_ptr1);
	temp = updateT( 0, copy_chosen_ptr1->parent );
	qA = qA * temp;


	selectOtherNodeP7( probs, nodes, reverse_copy_node_ptrs, reverse_copy_chosen_ptr2,
					   first_markers, theta, rho, locations );
	it = find( nodes.begin(), nodes.end(), reverse_copy_sib_ptr2);
	if ( it == nodes.end() ){
		/*cerr*/ Rcpp::Rcout << "ERROR: Reverse change is not possible in Major Haplotype update" << endl;
		throw Rcpp::exception("ERROR: Reverse change is not possible in Major Haplotype update"); //exit(1);
	}
	qA = qA * probs.at(int( it - nodes.begin() ));

	temp = updateT( 0, copy_chosen_ptr2->parent );
	qA = qA * temp;


	// Compute the values for the other updates but use the copy of the tree
	// because it relies on the old topology
	temp = updateUpTree( 0, copy_chosen_ptr1, copy_chosen_ptr2,
						copy_sib_ptr1, copy_sib_ptr2,
						copy_other_ptr1, copy_other_ptr2,
						rho, theta, first_markers, locations, ref_marker,
						allele_probs, haplo_probs );
	qA = qA * temp;


	deleteTree( reverse_copy_node_ptrs );

	sortNodePtrs(node_ptrs);

	// 4. Compute Pr(A~|G)
	int min_index = min( min(other_ptr1->id, sib_ptr1->id),
						 min(other_ptr2->id, sib_ptr2->id));

	if ( min_index < n ){
		min_index = n;
	}

	if (move_head == 1){
		max_ptr = node_ptrs[node_ptrs.size()-1];
	}

	probAnew = probP7_major ( chosen_ptr1,  sib_ptr1, other_ptr1,
							  chosen_ptr2,  sib_ptr2, other_ptr2,
		 					  node_ptrs, theta, rho, first_markers, ref_marker,
   							  locations,allele_probs, haplo_probs,
   							  min_index, max_ptr->id-1 );


	// 5. Compute Pr(A|G) using the copy of the tree that was made
	probA = probP7_major ( copy_chosen_ptr1, copy_other_ptr1, copy_sib_ptr1,
						   copy_chosen_ptr2, copy_other_ptr2, copy_sib_ptr2,
		 				   copy_node_ptrs, theta, rho, first_markers, ref_marker,
   						   locations,allele_probs, haplo_probs,
   						   min_index, max_ptr->id-1 );



	MHratio = min( exp(probAnew-probA)*exp(log(qA)-log(qAnew)), double( 1 ) );


	if ( myunif > MHratio ){

		deleteTree( node_ptrs[node_ptrs.size()-1] ) ;
		node_ptrs.clear();
		node_ptrs = copy_node_ptrs;
		accept = 0;

	} else {

		// Make the id of the first haplotype be the larger one
		if (node_ptrs[samples[chosen_index].seq1]->seq[samples[chosen_index].hetsites[0]]=='0'){
			swapNodes(samples[chosen_index], node_ptrs);
		}
		deleteTree( copy_node_ptrs[copy_node_ptrs.size()-1] );
		copy_node_ptrs.clear();
		updateflag=1;

	}

	return( accept );

}


void selectOtherNodeP7( vector<double>& probs, vector<treeNode*>& values,
					    vector<treeNode*>& node_ptrs, treeNode* chosen_node,
					    vector<int>& first_markers,double theta, double rho,
						   const vector<double>& locations, treeNode* other_node )
{

	unsigned int i=0;
	bool good_node=1;
	double total=0, epsilon=0.000001, lambda=0, t=0, prob;
	int j=0, l_index=0, r_index=0;

	probs.clear();
	values.clear();


	for (i=0; i<node_ptrs.size(); i++){

		// Determine if this node is eligible to be chosen.
		// Current node is the same as either of the two nodes chosen to
		// move. We can't connect a node to itself and we don't want to connect
		// it to the exact same spot
		if ( node_ptrs[i] == chosen_node ){
			good_node=0;
		}

		// This node is the parent of the chosen node and therefore cannot be
		// selected
		if ( node_ptrs[i] == chosen_node->parent ){
			good_node=0;
		}

		if ( other_node!=NULL ){
			if (node_ptrs[i] == other_node){
				good_node=0;
			}
		}


		if (good_node == 1){
			values.push_back(node_ptrs[i]);

			// Estimate the branch length separating these two nodes
			if (node_ptrs[i]->parent!=NULL){// Sampled from a uniform, use mean
				t = (node_ptrs[i]->parent->t + max(chosen_node->t,node_ptrs[i]->t))/2;
			} else {// Sampled from exp(1), use mean
				t = node_ptrs[i]->t+1;
			}
			lambda=(2*t-chosen_node->t-node_ptrs[i]->t)*theta/2;

			// Use the maximum amount of sequence that could be passed
			// to both nodes from their shared parent. Since chosen_node
			// is a tip node, this maximum is determined by its putative
			// sibling
			l_index=first_markers[0]-node_ptrs[i]->z_L+2;
			r_index=first_markers[1]+node_ptrs[i]->z_R-2;

			if ( (l_index>first_markers[0])&&(r_index<first_markers[1]) ){
				// If the sequence is not passed down to its ancestors at all due
				// to recombination events on either side, then there is no
				// sequence to compare, so let the probability of being chosen be
				// small (epsilon)
				probs.push_back(epsilon);

			} else {

				prob=1;
				for (j=l_index; j<=r_index; j++){
					if (node_ptrs[i]->seq[j] == chosen_node->seq[j]){
						prob = prob*(1/ALPHA*(1-exp(-lambda))+exp(-lambda));
					} else {
						prob = prob*1/ALPHA*(1-exp(-lambda));
					}
				}
				probs.push_back(prob);

			}
			total += probs[probs.size()-1];

		}

		good_node=1;

	}

	for (i=0; i<probs.size(); i++){
		probs[i] = probs[i]/total;

	}

}




double updateUpTree( bool sample, treeNode* chosen_1, treeNode* chosen_2,
					 treeNode* sib_1, treeNode* sib_2,
	  				 treeNode* other_1, treeNode* other_2,
	  				 const double rho, const double theta,
					 vector<int>& first_markers,
	  				 vector<double>& locations, int ref_location,
					 vector<vector<double> >& allele_probs,
	  				 vector<vector<double> >& hap_probs )
{

	unsigned int i=0, j=0, min=0;
	bool same=0;
	double prob=1/*, temp=0*/;


	
	// Create a vector containing the parental nodes. Remove any
	// duplicates if there are any
	treeNode* nodesarray[] = { chosen_1, chosen_2, other_1, other_2,
					           sib_1, sib_2 };
	vector<treeNode*> nodes(0);

	nodes.push_back(nodesarray[0]);

			  
	for (i=1; i<6; i++){
		for (j=0; j<nodes.size(); j++){
			if ( nodes[j] == nodesarray[i] ){
				same = 1;
			}
		}
		if (same != 1){
			if (nodesarray[i]!=NULL){
				nodes.push_back(nodesarray[i]);
			}
		}
		same=0;
	}

	vector<int> flags(nodes.size(),1);

	
	// Loop to update the parental information
	// NOTE: Shouldn't this be nodes.size()>=1???
	while ( (nodes[0]!=NULL) || (nodes.size()>1) ){

		// Find the element of nodes with minimum time
		min=0;
		for (i=1; i<nodes.size();i++){
			if ( nodes[i]->t < nodes[min]->t ){
				min=i;
			}
		}

		// Update the r and s value for this node or just find the probability
		// of the current values depending on whether we are sampling or not

		
			if ( nodes[min]->parent != NULL ){
				// We are updating using the same distribution from Pr(R|.) so
				// there is a cancellation in ratio.
				/*temp = */updateR( sample, nodes[min], rho, first_markers,
									   locations, ref_location, 4 );
				//prob = prob * temp;
				
			} else {
				nodes[min]->r_L=0;
				nodes[min]->r_R=0;
			}
			if ( flags[min] == 1 ){
				flags[min] = 0;
			} else {
				prob = prob * updateS_P7( sample, nodes[min], first_markers, theta,
						allele_probs, hap_probs );
			}


		// If the minimum's nodes parent is a duplicate, delete this node
		// from the vector of nodes
		same=0;
		i=0;
		if (nodes.size()>1){
			while ( (same==0)&&i<nodes.size() ){
				if ( nodes[i]==nodes[min]->parent ){
					same=1;

					if (flags[i]==1) {
						flags[i]=0;
					}
					nodes.erase(nodes.begin()+min);
					flags.erase(flags.begin()+min);
				}
				i++;
			}
		}

		// Move to the parent of the minimum provided that it wasn't already
		// in the vector of nodes.
		if ( same==0 ){
			nodes[min] = nodes[min]->parent;
		}			
	}

	if (nodes.size()>1){
		/*cout*/ Rcpp::Rcout << "ERROR in P7 - updateUpTree" << endl;
	}
	
	return(prob);
}

	


// This function is used to make the topology change for the sequence update
void topologyChangeP7( treeNode* node_ptr,
					   treeNode* other_ptr )
{

	// Determine the sibling of the node to be moved.
	// NOTE: Do this in here because two topology changes are done
	//       consecutively and the sib may have changed
	treeNode* sib_ptr = NULL;
	
	if (node_ptr->parent->child1->id == node_ptr->id){
		sib_ptr = node_ptr->parent->child2;
	} else {
		sib_ptr = node_ptr->parent->child1;
	}
	
	// Change the parent of the sibling to its current grandparent, and
	// change the grandparents child to the sibling
	sib_ptr->parent = sib_ptr->parent->parent;
	if (node_ptr->parent->parent!=NULL){ // In case the parent is the head node
		if (node_ptr->parent->parent->child1 == node_ptr->parent){
			sib_ptr->parent->child1 = sib_ptr;
		} else {
			sib_ptr->parent->child2 = sib_ptr;
		}
	}
	
	
	// Change the other node's parent to the parent of the node that is moving
	if (other_ptr->parent != NULL){
		if (other_ptr->parent->child1 == other_ptr){
			other_ptr->parent->child1 = node_ptr->parent;
		} else {
			other_ptr->parent->child2 = node_ptr->parent;
		}
	}

	
	// Make p's parent c's current parent
	node_ptr->parent->parent = other_ptr->parent;

	// Change c's parent to be p
	other_ptr->parent = node_ptr->parent;

	// Make c one of p's children
	if (node_ptr->parent->child1 == node_ptr){
		node_ptr->parent->child2 = other_ptr;
	} else {
		node_ptr->parent->child1 = other_ptr;
	}

	
	// Update the branch length for sib_ptr because it has a new parent
	// The time of the new parent has yet to be drawn so no need to update
	// the other nodes yet
	// NOTE: In order to use this function to reverse the topology change,
	//       all times should be adjusted. Note though that otherwise the
	//       tree will be nonsensical because the parent time has not been
	//       drawn to be less than its parent.
	if (sib_ptr->parent!=NULL){
		sib_ptr->b = sib_ptr->parent->t-sib_ptr->t;
	}
	if (node_ptr->parent->parent!=NULL){
		node_ptr->parent->b = node_ptr->parent->parent->t-node_ptr->parent->t;
	}
	node_ptr->b = node_ptr->parent->t-node_ptr->t;
	other_ptr->b = other_ptr->parent->t-other_ptr->t;

}


// This function computes a term proportional to the  log of the likelikood
// of the augmented data for the 7th proposal chain
double probP7_major( treeNode* node_ptr1, treeNode* other_ptr1, treeNode* sib_ptr1,
			   treeNode* node_ptr2, treeNode* other_ptr2, treeNode* sib_ptr2,
			   vector<treeNode*> node_ptrs, double theta, double rho,
	  		   vector<int>& first_markers, int ref_marker,
			   vector<double>& locations, vector<vector<double> >& allele_probs,
	   		   vector<vector<double> >& hap_probs,
	   		   int min_index, int max_index )
{
	double prob = log(double(1)), temp=0;
	unsigned int i=0, j=0, minimum=0;
	bool same=0;
	
	// Compute the prob(omega) terms.
	// NOTE: probT is already logged because some of the terms can get big
	temp =  probT( node_ptrs, min_index, max_index );
	prob = prob + temp;



	// Compute the prob(S) and prob(R) terms for the 6 lineages up the tree,
	// this is done similarly to the updates to ensure that the terms are not
	// duplicated
	treeNode* nodesarray[] = { node_ptr1, node_ptr2, sib_ptr1, sib_ptr2,
							   other_ptr1, other_ptr2 };
	vector<treeNode*> nodes(0);


	nodes.push_back(nodesarray[0]);

	for (i=1; i<6; i++){
		for (j=0; j<nodes.size(); j++){
			if ( nodes[j] == nodesarray[i] ){
				same = 1;
			}
		}
		if (same != 1){
			if (nodesarray[i]!=NULL){
				nodes.push_back(nodesarray[i]);
			}
			
		}
		same=0;
	}


	// Loop to compute the probabilities of the r and s values
	while ( (nodes[0]!=NULL) || (nodes.size()>1) ){

		// Find the element of nodes with minimum time
		minimum=0;
		for (i=1;i<nodes.size();i++){
			if ( nodes[i]->t < nodes[minimum]->t ){
				minimum=i;
			}
		}

		if ( nodes[minimum]->parent != NULL ){
			// We are updating using the same distribution from Pr(R|.) so
			// there is a cancellation in ratio.
			temp = probRi( nodes[minimum], rho, first_markers, locations,
					ref_marker );
			//prob = prob + log(temp);

		}
		temp = probSi( nodes[minimum], theta, first_markers, allele_probs,
					   hap_probs );

		prob = prob + log( temp );


		// If the minimum's nodes parent is a duplicate, delete this node
		// from the vector of nodes
		same=0;
		i=0;
		if (nodes.size()>1){
			while ( (same==0)&&i<nodes.size() ){
				if ( nodes[i]==nodes[minimum]->parent ){
					same=1;
					nodes.erase(nodes.begin()+minimum);
				}
				i++;
			}
		}

		// Move to the parent of the minimum provided that it wasn't already
		// in the vector of nodes.
		if ( same==0 ){
			nodes[minimum] = nodes[minimum]->parent;
		}
	}

	if (nodes.size()>1){
		/*cout*/ Rcpp::Rcout << "ERROR - In probP7" << endl;
	}

	return( prob );
	

}



// This  function is used for the proposal portion associated with a new
// sequence. For each locus it determines the PMF. Then either a new value
// for that locus is sampled based on the PMF and the probability is returned
// or the probability of the actual value, with respect to the PMF is returned.
double updateS_P7(bool sample, treeNode* node_ptr,
				  vector<int>& first_markers, double theta,
	  			  vector<vector<double> >& allele_probs,
				  vector<vector<double> >& hap_probs)
{
	int locus=0, locus_min=0, s_index=0, locus_max=0;
	int j=0;
	double prob=1;
	char temp=' ';
	vector<double> probs;

	// Update the left side
	locus_min = first_markers[0]-node_ptr->z_L+1;
	for ( locus = first_markers[0]; locus > locus_min; locus-- ){

		getSijPMF_P7("left", node_ptr, probs, first_markers, locus, theta,
					 allele_probs, hap_probs );


		// We are either sampling a new value for this locus based on Q()
		// (sample=1) and must also return the prob. Or we are computing the
		// probability of the reverse proposal (ie finding Q() for the allele
		// that we have at this locus
		if ( sample == 1 ){
			s_index = samplePMF( probs );
			changeSij( node_ptr, itoc( s_index ), locus );
			prob = prob * probs[s_index];
			if (prob==0){
				/*cout*/ Rcpp::Rcout << "Prob = 0; Left; sample" << endl;
			}
		} else {
			temp = node_ptr->seq[locus];
			if (temp == '0'){
				prob = prob * probs[0];
			} else if (temp == '1'){
				prob = prob * probs[1];
			} else { /*cout*/ Rcpp::Rcout << "Neither 0 nor 1" << endl; }
			if (prob==0){
				/*cout*/ Rcpp::Rcout << "Prob = 0; Left; Prob" << endl;
			}
		}

	}

	if ( sample == 1 ){
		for (j=locus; j>=0; j--){
			node_ptr->seq[j]='-';
		}
	}

	// Update the right side
	locus_max = first_markers[1]+node_ptr->z_R-1;
	for ( locus = first_markers[1]; locus < locus_max; locus++ ){

		getSijPMF_P7("right", node_ptr, probs, first_markers, locus, theta,
				  allele_probs, hap_probs );

		// We are either sampling a new value for this locus based on Q()
		// (sample=1) and must return the prob. Or we are computing the
		// probability of the reverse proposal (ie finding Q() for the allele
		// that we have at this locus
		if ( sample == 1 ){
			s_index = samplePMF( probs );
			changeSij( node_ptr, itoc( s_index ), locus );
			prob = prob * probs[s_index];
			if (prob==0){
				/*cout*/ Rcpp::Rcout << "Prob = 0; Right; sample" << endl;
			}
		} else {
			temp = node_ptr->seq[locus];
			if (temp == '0'){
				prob = prob * probs[0];
			} else if (temp == '1'){
				prob = prob * probs[1];
			} else { /*cout*/ Rcpp::Rcout << "Neither 0 nor 1" << endl; }
			if (prob==0){
				/*cout*/ Rcpp::Rcout << "Prob = 0; Right; prob" << endl;
			}
		}

	}

	if ( sample == 1 ){
		for (j=locus; j<(int)node_ptr->seq.size(); j++){
			node_ptr->seq[j]='-';
		}
	}

	return(prob);

}


// This finds the probability distribution for an s_j^i
void getSijPMF_P7(string side, treeNode* node_ptr, vector<double>& probs,
			   vector<int>& first_markers, int locus, double theta,
	  		   vector<vector<double> >& allele_probs,
   			   vector<vector<double> >& hap_probs)
{

	int sij=0;
	unsigned int i=0;
	string seq = node_ptr->seq;
	double probA=1, probC1=1, probC2=1, sum=0;
	int rC1=0, rC2=0;

	probs.clear();

	// Find the index (with respect to the sequence) for the recombination
	// variables for both of the children of node_ptr
	if ( side == "left" ){
		rC1 = first_markers[0]-node_ptr->child1->r_L+1;
		rC2 = first_markers[0]-node_ptr->child2->r_L+1;
	} else {
		rC1 = first_markers[1]+node_ptr->child1->r_R-1;
		rC2 = first_markers[1]+node_ptr->child2->r_R-1;
	}


	for (sij=0; sij<2; sij++){

		// Change the locus to either 0 or 1
		node_ptr->seq[locus]=itoc(sij);

		// Find probability associated with similarity to c1
		if ( ((side=="left")&&(locus>rC1)) || ((side=="right")&&(locus<rC1)) ){
			probC1 = probSij(side, node_ptr->child1, locus, theta, first_markers,
							 allele_probs, hap_probs);
		}

		// Find probability associated with similarity to c2
		if ( ((side=="left")&&(locus>rC2)) || ((side=="right")&&(locus<rC2)) ) {
			probC2=probSij(side, node_ptr->child2, locus, theta, first_markers,
						   allele_probs, hap_probs);
		}

		// The term corresponding to similarity with the parent. In this case
		// the parent will be modified so we will use information from this
		// node's sibling if they share this locus, or the haplotype model
		// if the don't
		probA=probSij_P7(side, node_ptr, locus, theta, first_markers,
				   allele_probs, hap_probs);

		probs.push_back(probA*probC1*probC2);
	}

	node_ptr->seq = seq;

	for (i=0; i<probs.size(); i++){
		sum = sum + probs[i];
	}
	for (i=0; i<probs.size(); i++){
		probs[i] = probs[i]/sum;
	}

}


// This computes the allele probability associated with a single locus. We pass
// the side and the value of j (1 to z-1). The probability comes from either
// the mutation model or the haplotype model, depending on j and z
double probSij_P7(string side, treeNode* node_ptr, int locus, double theta,
			   vector<int>& first_markers, vector<vector<double> >& allele_probs,
	 		   vector<vector<double> >& hap_probs)
{
	double probA=1;
	int prev_col=-1, prev_hcol=1;
	string prev_haplo(1,'0');
	int prev_row=0, cur_row=0;
	int ri=0, rs=0, r=0, type=1;
	treeNode* sib_ptr = NULL;
	double lambda = 0;

	// Determine sibling of this node if we aren't at the head node
	if (node_ptr->parent!=NULL){
		if (node_ptr->parent->child1 == node_ptr){
			sib_ptr = node_ptr->parent->child2;
		} else {
			sib_ptr = node_ptr->parent->child1;
		}

		lambda =  (node_ptr->b+sib_ptr->b)*theta/2;
	}


	// Values for the variables depends on whether we are on the left
	// or right side of the focal point
	if ( side == "left" ){

		ri = first_markers[0]-node_ptr->r_L+1;
		if ( (node_ptr->parent!=NULL) && (sib_ptr->r_L!=0) ){
			rs = first_markers[0]-sib_ptr->r_L+1;
			r = min(ri,rs);
		} else {
			r=ri;
		}


		// Determine whether this locus was inherited from this node's
		// parent (type=1) or not. If not, are we at the first locus to
		// recombine in (type=2) or not (type=3).
		if ( locus > r ){
			type = 1;
		} else {
			type = 2;
			prev_col = locus+1;
			prev_hcol = locus;
			if ( ((node_ptr->parent == NULL)&&(locus != first_markers[0]))||
					((node_ptr->parent != NULL)&&(locus < r)) ) {
				type = 3;
				prev_haplo = node_ptr->seq.substr(locus,2);
			}

		}


	} else {

		ri = first_markers[1]+node_ptr->r_R-1;
		if ( (node_ptr->parent!=NULL) && (sib_ptr->r_R!=0) ){
			rs = first_markers[1]+sib_ptr->r_R-1;
			r = min(ri,rs);
		} else {
			r=ri;
		}

		// Determine whether this locus was inherited to both node_ptr and its
		// sib from its parent (type=1). If not, are we at the first locus to
		// recombine in (type=2) or not (type=3).
		if ( locus < r ){
			type = 1;
		} else {
			type = 2;
			prev_col = locus-1;
			prev_hcol = locus-1;
			if ( ((node_ptr->parent == NULL)&&(locus != first_markers[1]))||
				((node_ptr->parent != NULL)&&(locus > r)) ){
				type = 3;
				prev_haplo = node_ptr->seq.substr(locus-1,2);
			}
		}

	}

	// This converts the two-locus haplotype sequence into a binary
	// variable. This avoids the use of if statements for "00", "01" etc.
	// Note that we must specify a default for next_haplo since it is in the
	// if statement. However, this variable is only used under the conditions
	// of the if, so if the value is not changed from the default, it doesn't
	// matter
	bitset<2> p_haplo (prev_haplo);


	// This sets the variables that index which rows of the allele/haplo
	// probability matrices. Converting in one step from a char element in
	// a string to a number ended up not working as expected.
	if (node_ptr->seq[prev_col] == '0'){
		prev_row = 0;
	} else { prev_row = 1; }
	if (node_ptr->seq[locus] == '0'){
		cur_row = 0;
	} else { cur_row = 1; }


	// Determine the prob of the given locus based on whether it is the first
	// locus after a recombination point (type 1) or not (type 2)
	if ( type == 1){ // Probability computed from mutation model
			if (node_ptr->seq[locus] == sib_ptr->seq[locus]){
				probA = 1/ALPHA*(1-exp(-lambda))+exp(-lambda);
			} else {
				probA = 1/ALPHA*(1-exp(-lambda));
			}
	} else if ( type == 2){ //Probability computed from haplo model
		probA=allele_probs[locus][cur_row];
	} else if (type == 3) {
		probA = hap_probs[prev_hcol][p_haplo.to_ulong()]/
				allele_probs[prev_col][prev_row];
	} else {
		/*cout*/ Rcpp::Rcout << "Problem: the type is neither 1, 2 or 3" << endl;
	}

	return(probA);

}



// The following function swaps information between the two nodes
// corresponding to a single individual
void swapNodes(genoNode chosen_node, vector<treeNode*>& node_ptrs){

	int tempint=0;
	treeNode* tempptr=NULL;

	// Swap the ids
	tempint=node_ptrs[chosen_node.seq1]->id;
	node_ptrs[chosen_node.seq1]->id=node_ptrs[chosen_node.seq2]->id;
	node_ptrs[chosen_node.seq2]->id=tempint;

	// Swap the pointers in node_ptrs
	tempptr=node_ptrs[chosen_node.seq1];
	node_ptrs[chosen_node.seq1]=node_ptrs[chosen_node.seq2];
	node_ptrs[chosen_node.seq2]=tempptr;
}


// This function is used if the dictionary rephaser is in use as
// well as the allele swap. If a new haplotype is added, this function
// compares it to the list of known haplotypes. If the new haplotype is
// not in the known list, this function will add the appropriate
// haplotype configurations.
void updateHaploList( const vector<string>& newseqs, vector<string>& haplo_list_vector,
		vector<genoNode>& samples, vector< vector<int> >& hapindexmat,
		vector<double>& P9scores )
{
  int id=0;
	unsigned int i=0, j=0, k=0/*, id=0*/;
	unsigned int numlist=haplo_list_vector.size();
	unsigned int numsamples=samples.size();
	unsigned int numtotalconfigs=hapindexmat.size();
	unsigned int numloci=samples[0].genotype.size();
	unsigned int numconfigs=0;
	bool goodhaplo=0, idfound=0, foundother=0;
	vector<int> found(2,0);
	vector<string> haplos(2);
	vector<int> addrow(2);
	string tempstr="";

	if ( (newseqs[0]=="")||(newseqs[1]=="") ){
		/*cout*/ Rcpp::Rcout << "A blank sequence" << endl;
	}

	// Check to see if BOTH of the two new haplotypes are in the list.
	// If both are, then nothing changes in the enumeration
	// Otherwise either or both of them will need to be added to the
	// list and the enumeration may need to be augmented.
	while ( (i<numlist)&&((found[1]==0)||(found[0]==0)) ){
		if (newseqs[0]==haplo_list_vector[i]){
			found[0]=1;
		} else if (newseqs[1]==haplo_list_vector[i]){
			found[1]=1;
		}
		i++;
	}

	// Add any haplotypes to the required lists
	for (i=0; i<2; i++){

		if (found[i]==0){

			// Add the haplotype to the haplotype list because it isn't
			// currently on there
			haplo_list_vector.push_back(newseqs[i]);


			// For each genotype in the set, compare the new haplotype to the
			// genotype. If the haplotype is not compatible, move on.
			// If it is compatible with the genotype, check its pair haplotype
			// (stored in haplos[1]). If the other one is already in the list
			// then this configuration will not need to be enumerated.
			haplos[0]=newseqs[i];


			for (j=0; j<numsamples; j++){

				goodhaplo=0;
				k=0;
				foundother=0;
				haplos[1]=haplos[0];
				idfound=0;

				// Check if this haplotype is compatible with the genotype
				// by comparing locus by locus. If the locus is compatible set
				// up the other haplotype
				while ( (goodhaplo==0)&&(k<numloci) ){

					if ( (samples[j].genotype[k]=='0')&&(newseqs[i][k]=='1') ){
						goodhaplo=1;
					} else if ( (samples[j].genotype[k]=='2')&&(newseqs[i][k]=='0') ){
						goodhaplo=1;
					} else {
						// Two haplotypes are initially equal, so only need to
						// modify the heterozygous loci
						if ( (samples[j].genotype[k]!='0')&&
								(samples[j].genotype[k]!='2') ){
							// Find which haplo should have the 0 or 1
							if (haplos[0][k]=='0'){
								haplos[1][k]='1';
							} else {
								haplos[1][k]='0';
							}
						}
					}
					k++;
				}


				// It is compatible. If the other haplotype of the pair is not
				// already in the list, add this configuration to the
				// enumeration related vectors.
				if ( goodhaplo==0 ){

					// Check whether the other haplotype is already in the list
					k=0;
					numlist=haplo_list_vector.size();
					while ( (k<numlist)&&(foundother==0)){
						if (haplos[1]==haplo_list_vector[k]){
							foundother=1;
						}
						k++;
					}

					if (foundother!=1){


						//first_het_site=samples[j].hetsites[0];
						//if (haplos[0][first_het_site]=='0'){
							// The first one should be the larger, so swap the two
						//	tempstr=haplos[0];
						//	haplos[0]=haplos[1];
						//	haplos[1]=tempstr;
						//}

						// Add it to the information for this sample
						numconfigs=samples[j].hapconfigs.size();
						numtotalconfigs=hapindexmat.size();
						samples[j].hapconfigs.push_back(haplos);
						id=samples[j].id;

						k=0;
						while ( (k<numtotalconfigs)&&(idfound==0) ){
							if (hapindexmat[k][0]==(id-1) ){
								idfound=1;
							}
							k++;
						}

						if ( (k==numtotalconfigs)&&(idfound==0) ){
							/*cerr*/ Rcpp::Rcout << "ERROR IN updateHaploList. Could not find ID" << endl;
							throw Rcpp::exception("ERROR IN updateHaploList. Could not find ID"); //exit(1);
						} else {
							addrow[0]=id-1;
							addrow[1]=numconfigs;
							hapindexmat.insert(hapindexmat.begin()+k+numconfigs-1,addrow);
							P9scores.insert(P9scores.begin()+k+numconfigs-1,0);
						}


					}
				}
			}


		}

	}


}
