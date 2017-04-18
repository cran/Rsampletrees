/*----------------------------------------------------------------------------

  readFiles.cpp - This contains the functions for reading in the data and 
                  options in the sequence data file and parameter file 
                  respectively. It also has function to print what was read
                  in and to summarize the information.
  
  Kelly Burkett; March 15, 2007

  Modification Notes:

  --------------------------------------------------------------------------*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <Rcpp.h>

#include "misc.h"
#include "nodefunctions.h"
#include "readFiles.h"
#include "treeBuild.h"
#include "rng.h"

using namespace std;

void readOptions(Options& myoptions, ifstream& myparamfile)
{
	
	string line="", option_name="", option_value="";
	int line_count=0/*, i=0, maxweights=5, flag=0*/;
//	double sum=0;
	
	while ( !myparamfile.eof() )
	{
		line_count += 1;
		myparamfile >> option_name;
		myparamfile >> option_value;
		
		// Take the options read in to variables option_name and option_value
		// and copy them to the appropriate field of myoptions
		if (option_name == "DataFile"){
			myoptions.datafile = option_value; }
		else if (option_name == "RunName"){
			myoptions.run_name = option_value;}
		else if (option_name == "LocationFile"){
			myoptions.location_file = option_value; }
		else if (option_name == "DetailedOutput"){
			myoptions.output = atoi(option_value.c_str()); }
		else if (option_name == "ChainLength"){
			myoptions.len_chain = atoi(option_value.c_str()); }
		else if (option_name == "Seed"){
			// ato- functions are not available for long long int's so
			// I use the stringstream trick.
			stoLL(option_value, myoptions.seed);
		}
		else if (option_name == "BurnIn"){
			myoptions.burn_in = atoi(option_value.c_str()); }
		else if (option_name == "FocalPoint"){
			myoptions.x = atoi(option_value.c_str()); }
		else if (option_name == "InitialTree"){
			myoptions.initial_tree = atoi(option_value.c_str());
			if ( (option_value!="0") && (option_value!="1") ){
				/*cerr*/ Rcpp::Rcout << "Invalid option for 'InitialTree'. Should be 0 (no) or 1 (yes)" << endl;
				throw Rcpp::exception("Invalid option for 'InitialTree'. Should be 0 (no) or 1 (yes)"); //exit(1);
			}
		}
		else if (option_name == "InitialTreeFile"){
			myoptions.initial_tree_file = option_value; }
		else if (option_name == "RandomTree"){
			myoptions.random = atoi(option_value.c_str()); }
/*		else if (option_name == "ST"){
			myoptions.st_rho = atoi(option_value.c_str()); }
		else if (option_name == "SA"){
			myoptions.st_SA = atoi(option_value.c_str()); }
		else if (option_name == "STFile"){
			myoptions.st_file = option_value; }
*/		else if (option_name == "InitialHaplos"){
			myoptions.initial_haplos = atoi(option_value.c_str());
			if ( (option_value!="0") && (option_value!="1") ){
				/*cerr*/ Rcpp::Rcout << "Invalid option for 'InitialHaplos'. Should be 0 (no) or 1 (yes)" << endl;
				throw Rcpp::exception("Invalid option for 'InitialHaplos'. Should be 0 (no) or 1 (yes)"); //exit(1);
			}
		}
		else if (option_name == "InitialHaploFile"){
			myoptions.initial_haplo_file = option_value; }
		else if (option_name == "DataType"){
			myoptions.datatype = *option_value.c_str();
			if ( (myoptions.datatype!='h') && (myoptions.datatype!='g') ){
				/*cerr*/ Rcpp::Rcout << "Invalid option for 'DataType'. Should be h (haplotype) or g (genotype)" << endl;
				throw Rcpp::exception("Invalid option for 'DataType'. Should be h (haplotype) or g (genotype)"); //exit(1);
			}
		}
		else if (option_name == "HaploFreqFile"){
			myoptions.haplo_freq_file = option_value; }
		else if (option_name == "InitialTheta"){
			myoptions.initial_theta = atof(option_value.c_str()); }
		else if (option_name == "MinTheta"){
			myoptions.min_theta = atof(option_value.c_str()); }
		else if (option_name == "MaxTheta"){
			myoptions.max_theta = atof(option_value.c_str()); }
		else if (option_name == "InitialRho"){
			myoptions.initial_rho = atof(option_value.c_str()); }
		else if (option_name == "ScaleRho"){
			myoptions.scale_rho = atof(option_value.c_str()); }
		else if (option_name == "ShapeRho"){
			myoptions.shape_rho = atof(option_value.c_str()); }
		else if (option_name == "Thinning"){
			myoptions.thinning = atoi(option_value.c_str()); }
		else if (option_name == "WeightFile"){
			myoptions.weight_file = option_value; }
		else if (option_name == "HaploList"){
			myoptions.haplo_list = atoi(option_value.c_str());
			if ( (option_value!="0") && (option_value!="1") ){
				/*cerr*/ Rcpp::Rcout << "Invalid option for 'HaploList'. Should be 0 (no) or 1 (yes)" << endl;
				throw Rcpp::exception("Invalid option for 'HaploList'. Should be 0 (no) or 1 (yes)"); //exit(1);
			}
		}
		else if (option_name == "HaploListFile"){
			myoptions.haplo_list_file = option_value; }
		else if (option_name == ""){}
		else {
			/*cerr*/ Rcpp::Rcout << "The options file has an unknown option on line ";
			/*cerr*/ Rcpp::Rcout << line_count << endl;
			throw Rcpp::exception("The options file has an unknown option on line"); //exit(1);
		}
		option_name="";	
		getline(myparamfile, line);
				
	}

}


// This function reads in the sequence data from the file specified in the 
// parameter file. It is read into the string vector seqs. Having a string element of 
// " " will cause problems later on, so it must be addressed by replacing 
// or aborting.
void readSeq(ifstream& myseqfile, vector<string>& seqs)
{
	string temp;
	int index=0;
	unsigned int seqlen=0;
	
	while ( !myseqfile.eof() )
	{
		getline (myseqfile,temp);

		if (!temp.empty()){ // To strip any blanks at the end of the file
			if ( index == 0 ) {
				seqlen = temp.length();
			}
		
			if ( temp.length() != seqlen ){
				/*cerr*/ Rcpp::Rcout << "Sequence on line " << index+1 << " does not have ";
				/*cerr*/ Rcpp::Rcout << seqlen << " elements\n";
				throw Rcpp::exception("Sequence does not have proper number of elements"); //exit(1);
			} else {
				seqs.push_back(temp);
			}
			index++;
		}
	}
	Rcpp::Rcout << "Sequence Information" << endl;
	Rcpp::Rcout << "\t\tNumber of sequences:" << seqs.size() << endl;
	Rcpp::Rcout << "\t\tSequence length:" << seqlen << endl << endl;
	
}

// This function is used to read the genotype data from a genotype file into
// a vector of genotype data nodes. 
void readGenos( ifstream& genofile, vector<genoNode>& genodat,
				vector<int>& hetIDs, vector<vector<int> >& hapindexmat,
				int dictflag, const vector<string>& haplolist)
{
	string temp="";
	string geno="";
	unsigned int i=0, numloci=0;
	string newvalue;
	genoNode newNode = {"",0,-1,-1};
	vector<string> haplos(2);
	vector<vector<string> > allhaplos;
	vector<int> aset(2,0);
	bool isMissing=0;

	while ( !genofile.eof() )
	{
		getline (genofile,temp);

		if ( !temp.empty() ){ // To strip any blanks at the end of the file

			genodat.push_back(newNode);
			genodat[genodat.size()-1].id = genodat.size();
			
			for (i=0; i<temp.size(); i++){
				newvalue = temp[i];
				if ( newvalue == " " ){
				} else if ( (newvalue=="0") || (newvalue=="1") ||
							(newvalue=="2") ){
					geno.push_back(temp[i]);
					if ( newvalue == "1"){
						(genodat[genodat.size()-1].hetsites).push_back((i+1)/2);
						
					}
				} else {
					/*cerr*/ Rcpp::Rcout << "ERROR - Invalid genotype on line "
						 << genodat.size();
					/*cerr*/ Rcpp::Rcout << " of genotype data file" << endl;
					throw Rcpp::exception("ERROR - Invalid genotype"); //exit(1);
				}
			}
			if ( genodat.size()==1 ){
				numloci = geno.size();
			}
			if ( geno.size() != numloci ){
				/*cerr*/ Rcpp::Rcout << "ERROR - Incorrect number of genotypes on line "
					 << genodat.size();
				/*cerr*/ Rcpp::Rcout << " of genotype data file" << endl;
				throw Rcpp::exception("ERROR - Incorrect number of genotypes"); //exit(1);
			}
			genodat[genodat.size()-1].genotype = geno;
			if ( genodat[genodat.size()-1].hetsites.size()>1 ){
				isMissing=1;
				hetIDs.push_back(genodat.size());
			}

			if (dictflag==1){
				if (haplolist.size()==0){
					findHaplos_all(haplos,allhaplos,
							genodat[genodat.size()-1].genotype,0,numloci,1);
				} else {
					if ( numloci!=haplolist[0].size() ){
						/*cerr*/ Rcpp::Rcout << "ERROR - Mismatch between number of loci in"
							 << "genotype file and in list of likely "
							 << "haplotypes" << endl;
						throw Rcpp::exception("ERROR - Mismatch between number of loci in genotype file and in list of likely haplotypes"); //exit(1);
					}
					findHaplos_approx(isMissing, allhaplos,
							genodat[genodat.size()-1],haplolist);
				}

				genodat[genodat.size()-1].hapconfigs=allhaplos;
				for (i=0; i<allhaplos.size(); i++){
					aset[0]=genodat[genodat.size()-1].id-1;
					aset[1]=i;
					hapindexmat.push_back(aset);
				}
				allhaplos.clear();
				haplos[0]="";
				haplos[1]="";
			}
			geno="";
			isMissing=0;
			
		}
		
	}
	Rcpp::Rcout << genodat.size() << " rows read in from the genotype data file"
		 << endl;


}



// This function reads the marker locations from the distance file into a vector
// of type double. In the future, the type may have to change depending on
// how distances will be read.
void readDist(ifstream& mydistfile, vector<double>& dists){
	
	unsigned int i=0;

	Rcpp::Rcout << "Marker locations:" << endl;
	Rcpp::Rcout << "\t\t";
	while ( (!mydistfile.eof()) && (i!=dists.size()) )
	{
		mydistfile >> dists[i];
		Rcpp::Rcout << dists[i] << " ";

		if ( (i>0)&&(dists[i] < dists[i-1]) ){
			// Locations were not in increasing order so fail
			/*cerr*/ Rcpp::Rcout << "Location " << i << " is less than location " << i-1 ;
			/*cerr*/ Rcpp::Rcout << ". Check file." << endl;;
			throw Rcpp::exception("Location problem. Check file."); //exit(1);
		}
		i++;
		
	}
	Rcpp::Rcout << endl;

	if ( i < dists.size() ){
		// Too few locations read in so fail
		/*cerr*/ Rcpp::Rcout << "The number of locations read in from the input file ";
		/*cerr*/ Rcpp::Rcout << "was less than the number of markers read in from the";
		/*cerr*/ Rcpp::Rcout << " the sequence file" << endl;
		throw Rcpp::exception("Number of locations read less than number of markers read."); //exit(1);
	}

	if ( !mydistfile.eof() ){
		// Too many locations could be read in so input was truncated.
		// This can occur both if there are too many locations or if there
		// are hidden blanks or end-of-line characters at the end of the file.
		// Only a warning is given.
		Rcpp::Rcout << endl << "WARNING: Only the first " << dists.size() << " locations ";
		Rcpp::Rcout << "read in from the location file." << endl;
		Rcpp::Rcout << "The location file may contain more locations";
		Rcpp::Rcout << " than expected, blank spaces or hidden characters";
		Rcpp::Rcout << endl << endl;
	}
}



// This function reads in a tree in Newick format, and creates a tree structure
// out of it. That tree can then be used to have the partitions found, or to
// find the average pairwise distance (once I can read in times that is)
treeNode* readTree(ifstream& treefile, int id)
{
	treeNode *head_node = createNode();
	treeNode *current_node ;
	treeNode *parent_node = head_node;
	head_node->id = id;
	id--;
	
	string tree, sub="", temp;

	// Read in the tree. It should be one line in the input file
	if ( !treefile.eof() ){
		getline ( treefile, tree );
	}

	// Parse the string that was read in from the file. 
	unsigned int i=1;
	current_node = NULL;
	while  (i<tree.size()) {

		temp = tree[i];
		
		if ( temp == "("  ){
			
			// Create a new internal node and set its parent (if applicable)
			// and set it as a child of its parent. 
			current_node = createNode();
			current_node->parent = parent_node;
			current_node->id = id;
			id--;

			if (parent_node->child1 == NULL){
				parent_node->child1 = current_node;
			} else {
				parent_node->child2 = current_node;
			}

			// Make the newly created node the new parent and move over 1
			parent_node = current_node;
			i++;
			
		} else if ( temp == ")" ){
			
			// Move up a level in the tree
			current_node = current_node->parent;
			parent_node = current_node->parent;
			i++;
			
		} else if ( temp == "," ){

			// Move over 1 since another node follows
			i++;
			
		} else if (temp == ":" )  {

			// Get the time of the most recent node
			i++;
			do {
				sub.push_back(tree[i]);
				i++;
				temp = tree[i];
			} while ( isdigit( tree[i] ) || (temp == "."));

			current_node->b = atof ( sub.c_str() );
			sub.clear();
			
		} else if (temp == ";" ) {
			// Move over 1 
			i++;
		} else {

			// Create a new terminal node, and set the parents/children
			// of this node and its parent
			current_node = createNode();
		
			current_node->parent = parent_node;
			if (parent_node->child1 == NULL){
				parent_node->child1 = current_node;
			} else {
				parent_node->child2 = current_node;
			}
			
			// Find the id of this node 		
			do {
				sub.push_back(tree[i]);
				i++;
			} while ( isdigit( tree[i]) );
			
			current_node->id = atoi ( sub.c_str() );
			sub.clear();
			
		}
		
	}
	getTimes( head_node );
	return( head_node );
}


// This function gets the times since the present from the branch length
// times that are available
void getTimes ( treeNode* node_ptr ){

	
	if ((node_ptr->child1->child1 == NULL)&&(node_ptr->child2->child1 == NULL)){
		
		// Left and Right grandchildren of current node are NULL (at bottom)
		node_ptr->t = node_ptr->child1->b;		

	} else if ((node_ptr->child1->child1 != NULL)&&(node_ptr->child2->child1 == NULL)){
		
		// Only Right grandchild is NULL. So at bottom on right side of tree, but
		// not left, so here "b" for the right child is the time of this node
		getTimes ( node_ptr->child1 );
		node_ptr->t = node_ptr->child2->b;

		
	} else if ((node_ptr->child1->child1 == NULL)&&(node_ptr->child2->child1 != NULL)){
		
		// Only Left grandchild is NULL. So at bottom on left side of tree, but
		// not right, so here "b" for the left child is the time of this node
		getTimes ( node_ptr->child2 );
		node_ptr->t = node_ptr->child1->b;
		
	} else {
		
		// Subtrees on both left and right hand side to print
		getTimes ( node_ptr->child1 );
		getTimes ( node_ptr->child2 );
		node_ptr->t = node_ptr->child1->b + node_ptr->child1->t;
		
	}

}




// This function summarizes the sequence data in a vector of type SeqData.
// Each element of SeqData contains a sequence, the number of individuals
// that have this sequence, and the ID (in order they were read in) of 
// the individuals having this sequence.
vector<SeqData> uniqueSeq(const vector<string> seqs){
	
	vector<SeqData> myseqinfo(1);
	SeqData new_info;
	unsigned int i, j, current_size;
	bool found=0;
	
	// Set up the first element of the vector myseqinfo
	myseqinfo[0].seq_type = seqs[0];
	myseqinfo[0].num_of_type = 1;
	myseqinfo[0].seq_ID[0] = 1;
	
	for (i=1; i<seqs.size(); i++){ // over the list of sequences in seqs[]
		
		current_size = myseqinfo.size();
		found = 0; // indicator of whether the sequence has been found in the list
		j = 0;     // which element of the vector we are at in our search
		
		while (!found){
			
			if (seqs[i] == myseqinfo[j].seq_type){ 
				
				// found the sequence in myseqinfo
				myseqinfo[j].num_of_type++;
				myseqinfo[j].seq_ID.push_back(i+1);
				found = 1;
				
			} else if ( j == (current_size-1) ){ 
				
				// sequence isn't in myseqinfo, so add a new element (of type
				// SeqData to the vector myseqinfo
				myseqinfo.push_back(new_info);
				myseqinfo[current_size].seq_type = seqs[i];
				myseqinfo[current_size].num_of_type = 1;
				myseqinfo[current_size].seq_ID[0] = (i+1);
				found = 1;
				
			}
			j++;
				
		}
			
	} 
	
	return(myseqinfo);
}


// This function mirrors back what the options were to screen so that 
// we can ensure that all values were input correctly.
void printOptions(Options& myoptions, ifstream& myparamfile)
{

	Rcpp::Rcout << "Data files" << endl;
	Rcpp::Rcout << "\t\tGenotype/Sequence data: " << myoptions.datafile << endl;
	Rcpp::Rcout << "\t\tMarker location data: " << myoptions.location_file << endl;
	if (myoptions.initial_tree == 1){
		Rcpp::Rcout << "\t\tInitial tree data:"<< myoptions.initial_tree_file << endl;
	} else {
		Rcpp::Rcout <<"\t\tNo initial tree data file specified" << endl;
		Rcpp::Rcout << "\t\t\t\tInitial data will be generated" << endl;
		if (myoptions.random == 0){
			Rcpp::Rcout << "\t\t\t\tTree generated via UPGMA" << endl;
		} else {
			Rcpp::Rcout << "\t\t\t\tTree generated by randomly connecting nodes" << endl;
		}
	}
	if (myoptions.datatype == 'g'){
		if (myoptions.initial_haplos == 1){
			Rcpp::Rcout << "\t\tInitial haplotype configuration data:"
				 << myoptions.initial_haplo_file << endl;
		} else {
			Rcpp::Rcout <<"\t\tNo initial haplotype configuration data file specified" << endl;
			Rcpp::Rcout << "\t\t\t\tInitial haplotype configuration data will be generated" << endl;
		}
		Rcpp::Rcout << "\t\tHaplotype frequency File: " << myoptions.haplo_freq_file << endl;
	}
	Rcpp::Rcout << endl;

	Rcpp::Rcout << "User Options" << endl;
	Rcpp::Rcout << "\t\tDetailed output to screen (0 no; 1 yes): " << myoptions.output << endl;
	Rcpp::Rcout << "\t\tGenotype or Haplotype (h haplotype; g genotype): "
		 << myoptions.datatype << endl;
	Rcpp::Rcout << "\t\tNumber of MCMC samples: " << myoptions.len_chain << endl;
	Rcpp::Rcout << "\t\tBurn-in samples: " << myoptions.burn_in << endl;
	Rcpp::Rcout << "\t\tThining interval: " << myoptions.thinning << endl;
	Rcpp::Rcout << "\t\tFocal point location: " << myoptions.x << endl;

	if (myoptions.initial_theta == -1){
		Rcpp::Rcout << "\t\tNo initial theta specified" << endl;
		Rcpp::Rcout << "\t\t\t\t Default theta to be used" << endl;
	}
	if (myoptions.min_theta == -1){
		Rcpp::Rcout << "\t\tNo minimum theta for uniform distribution specified" << endl;
		Rcpp::Rcout << "\t\t\t\tDefault minimum theta to be used" << endl;
	}
	if (myoptions.max_theta == -1){
		Rcpp::Rcout << "\t\tNo maximum theta for uniform distribution specified" << endl;
		Rcpp::Rcout << "\t\t\t\tDefault maximum theta to be used" << endl;}
	if (myoptions.initial_rho == -1){
		Rcpp::Rcout << "\t\tNo initial rho specified" << endl;
		Rcpp::Rcout << "\t\t\t\tDefault rho to be used" << endl;}
	if (myoptions.scale_rho == -1){
		Rcpp::Rcout << "\t\tNo scale parameter for the gamma prior for rho" << endl;
		Rcpp::Rcout << "\t\t\t\tDefault to be used" << endl;}
	if (myoptions.shape_rho == -1){
		Rcpp::Rcout << "\t\tNo shape parameter for the gamma prior for rho" << endl;
		Rcpp::Rcout << "\t\t\t\tDefault to be used" << endl;}
	Rcpp::Rcout << endl;


}


// This function is used to display the list of sequences to ensure that 
// no input errors were made.
void printSeq(const vector<string>& seqs)
{
	unsigned int i;

	Rcpp::Rcout << "Sequence Data:" << endl;
	for (i=0; i<seqs.size(); i++){
		Rcpp::Rcout << "\t\t" << seqs[i] <<endl;
	}
	Rcpp::Rcout << endl;
	
}


// This function reads in the haplotype frequencies from a file specified by the user
void readHaps( ifstream& hapfile, vector<vector<double> >& hap_freqs, unsigned int numcols )
{
	int i=0, j=0;
	unsigned int counter = 0;
	double sum=0;
	vector<double> row(4,-1);

	while ( (!hapfile.eof()) && (counter < numcols*4) ){
 
		if ( j==0 ){
			hap_freqs.push_back(row);
		}
		hapfile >> hap_freqs[i][j];
		sum += hap_freqs[i][j];
		
		j++;
		if ( j==4 ){ //Verify contents of this row and start a new one
			if ( abs(sum-1) >=0.00001 ){
				/*cerr*/ Rcpp::Rcout << "ERROR - Haplotype proportions do not"
						<< " sum to 1 in line " << i+1
						<< " of frequency file" << endl;
				throw Rcpp::exception("Haplotype proportions do not sum to 1"); //exit(1);
			}
			i++;
			j=0;
			sum=0;
			
		}
		counter++;
		
	}

	if ( (hap_freqs.size()<numcols) || hap_freqs[numcols-1][3] == -1 ){
		/*cerr*/ Rcpp::Rcout << "ERROR - Too few values in haplotype"
				<<	" probability file" << endl;
		throw Rcpp::exception("Too few values in haplotype"); //exit(1);
	}

	if ( !hapfile.eof() ){
		Rcpp::Rcout    << endl << "WARNING - Only the first " << counter
				<< " values read in from the haplotype file."
				<< " File may contain too many values or "
				<< "hidden characters" << endl << endl;
	}

        
}


// This function sets up an initial set of sequences based on the genotype
// data and the haplotype frequencies
void initialHaplos_freqs( vector<string>& seqs, vector<genoNode>& data,
				    const vector<vector<double> >& hap_freqs )
{

	unsigned int i=0, j=0;
	int allele=0;
	string seq1="";
	string seq2="";
	string geno="";
	int seqflag=0;

	for (i=0; i<data.size(); i++){

		geno = data[i].genotype;

		for (j=0; j<geno.size(); j++){

			if ( geno[j] == '0') {
				seq1.push_back('0');
				seq2.push_back('0');
			} else if ( geno[j] == '2' ) {
				seq1.push_back('1');
				seq2.push_back('1');
			} else {
				if ( (j==0) || (geno[j-1]=='0') || (geno[j-1]=='2') ){
					allele = generateDiscreteUnif(0,1);
					seq1.push_back(itoc(allele));
					seq2.push_back(itoc(1-allele));
				} else {
					if (seq1[seq1.size()-1] =='1'){
						seqflag=1;
					}
					if ( (hap_freqs[j-1][0]*hap_freqs[j-1][3])>=(hap_freqs[j-1][1]*hap_freqs[j-1][2]) ){
						if (seqflag==0){
							seq1.push_back('0');
							seq2.push_back('1');
						} else {
							seq1.push_back('1');
							seq2.push_back('0');
						}
					} else {
						if (seqflag==0){
							seq1.push_back('1');
							seq2.push_back('0');
						} else {
							seq1.push_back('0');
							seq2.push_back('1');
						}
					}
				}
				
			}
		}

		if (data[i].hetsites.size()==0){
			seqs.push_back(seq1);
			seqs.push_back(seq2);
		} else {
			if (seq1[data[i].hetsites[0]]=='1'){

				// The first sequence has a higher decimal number
				seqs.push_back(seq1);
				seqs.push_back(seq2);

			} else {

				// The second sequence has a higher decimal number
				seqs.push_back(seq2);
				seqs.push_back(seq1);

			}
		}
		data[i].seq1 = seqs.size()-2;
		data[i].seq2 = seqs.size()-1;
		seq1="";
		seq2="";

	}

}


// This function sets up an initial set of sequences based on the genotype
// data and the haplotype frequencies
void initialHaplos_enum( vector<string>& seqs, vector<genoNode>& data )
{

	unsigned int i=0, num=0, sample=0;
	string seq1="";
	string seq2="";
	string geno="";

	for (i=0; i<data.size(); i++){

		// Sample a configuration from the enumeration
		num=data[i].hapconfigs.size();
		sample = generateDiscreteUnif(0,num-1);
		seq1=data[i].hapconfigs[sample][0];
		seq2=data[i].hapconfigs[sample][1];

		// Add the sampled configuration to the sequence list. Order the
		// higher "number" first.
		if (data[i].hetsites.size()==0){
			seqs.push_back(seq1);
			seqs.push_back(seq2);
		} else {
			if (seq1[data[i].hetsites[0]]=='1'){
				// The first sequence has a higher decimal number
				seqs.push_back(seq1);
				seqs.push_back(seq2);

			} else {
				// The second sequence has a higher decimal number
				seqs.push_back(seq2);
				seqs.push_back(seq1);

			}
		}

		// Store the label for the two haplotypes corresponding to id 'i'
		data[i].seq1 = seqs.size()-2;
		data[i].seq2 = seqs.size()-1;
	}

}

void readst( vector<double>& constants, vector<double>& ranges,
			 ifstream& stfile, Options myoption)
{
	double value=0, shape=0;
	int line=0/*, i=0*/;

	while ( stfile >> value ){

		line++;
		constants.push_back(value);

		stfile >> shape;

		if (shape<0){
			/*cerr*/ Rcpp::Rcout << "Problem with the shape parameter given for rho "
				 << "on line " << line << " of the ST file" << endl;
			throw Rcpp::exception("Problem with the shape parameter given for rho"); //exit(1);
		}

		ranges.push_back(shape);

	}

}



void readWeights( vector<double>& weights, vector<vector<int> >& updates,
		          string wfile, char type, /*bool st,*/ int& dictflag)
{
	double defweights[]={0.05,0.05,0.5,0.15,0.15,0.1,0};
	ifstream weightfile;
	vector < int > myints;
	double mydbl;
	int myint;
	string line;
	int i=0, j=0;
	double sum=0;
	//int STflag=0;

	if ( wfile=="" ){
		// Use default weights
		weights.insert ( weights.begin(), defweights, defweights+7 );
		sum=1;

		if ( type=='g' ){// genotype data
			weights[6]=0.1;
			weights[4]=0.1;
			weights[3]=0.1;
		}
/*
		if ( st==1 ){// do ST so reweight the weights
			weights.push_back(0.1);
			sum=0.1;
			for (i=0; i<(weights.size()-1); i++){
				weights[i] = weights[i]*(1-0.1);
				sum += weights[i];
			}
		}
*/
		for (i=1; i<=(int)weights.size(); i++){
			myints.push_back(i);
			updates.push_back(myints);
			myints.clear();
		}


	} else {
		// Weights provided in a file by the user
		weightfile.open( wfile.c_str(), ios::in );
		if ( weightfile.fail() ){
			/*cerr*/ Rcpp::Rcout << "The weights file " << wfile << " cannot be opened."
				 << endl;
			throw Rcpp::exception("The weights file cannot be opened."); //exit ( 1 );
		}

		while ( getline(weightfile, line) ){

			istringstream iss(line);

			// The first number is the weight associated with the update
			iss >> mydbl;
			sum+=mydbl;
			weights.push_back(mydbl);
			mydbl=0;

			while (iss >> myint){
				myints.push_back(myint);
				if (myint==7){
					// A dictionary rephaser update is being done
					dictflag=1;
				}
/*				if ((myint==8)&&(st==0)){
					Rcpp::Rcout << "ERROR: ST option turned off but weights file "
							"gives weight to update 8\n" << endl;
					throw Rcpp::exception("ST option turned off but weights file gives weight to update 8"); //exit(1);
				}
				if ((myint==2)&&(st==1)){
					Rcpp::Rcout << "WARNING: Rho updates (2) given positive probability and ST is on. "
							"No Rho updates will be done except with ST." << endl;
				} */
//				if (myint==8){ STflag=1; }
			}
			if (myints.size()==0){
				/*cerr*/ Rcpp::Rcout << "ERROR: Problem in weights file. Need at least two "
						"values per line\n" << endl;
				throw Rcpp::exception("Problem in weights file"); //exit(1);
			}

			updates.push_back(myints);
			myint=0;
			myints.clear();
		}

		weightfile.close();
	}

	Rcpp::Rcout << "\t\tWeights and corresponding update types" << endl;
	for (i=0; i<(int)weights.size(); i++){
		Rcpp::Rcout << "\t\t\t\t" << weights[i] << " (";
		for (j=0; j<(int)updates[i].size(); j++){
			Rcpp::Rcout << updates[i][j];
			if (j!=((int)updates[i].size()-1)){
				Rcpp::Rcout << ",";
			}
		}
		Rcpp::Rcout << ")" << endl;
	} Rcpp::Rcout << endl;
/*
	if ( (st==1) && (STflag==0) ){
		Rcpp::Rcout << "ERROR: ST option turned on but weights file gives 0 weight to"
				" update 8\n" << endl;
		throw Rcpp::exception("ST option turned on but weights file gives 0 weight to update 8"); //exit(1);
	}
*/
	if (abs(sum-1)>0.0000001){
		/*cerr*/ Rcpp::Rcout << "ERROR: weights do not sum to 1:" << sum << endl;
		throw Rcpp::exception("ERROR: weights do not sum to 1"); //exit(1);
	}
}


// If the dictionary rephaser is being done the list of haplotype
// configurations must be put together. If a file of likely
// haplotypes is not provided, all configurations are enumerated
// by the function findHaplos_all(). If the list is provided, an
// approximation consisting of only combinations of haplotypes in the
// list is used. This approximation is found with the function
// findHaplos_approx()
void findHaplos_all(vector<string>& haplos, vector<vector<string> >& allhaps,
		string genotype, int locus, int numloci, bool firsthet)
{

	string haplo0, haplo1;

	if (locus<numloci){

		if (genotype[locus]=='0'){
			haplos[0].push_back('0');
			haplos[1].push_back('0');
			findHaplos_all(haplos,allhaps,genotype,locus+1,numloci,firsthet);
		} else if (genotype[locus]=='2'){
			haplos[0].push_back('1');
			haplos[1].push_back('1');
			findHaplos_all(haplos,allhaps,genotype,locus+1,numloci,firsthet);
		} else {
			if (firsthet==1){ // In the enumeration, the first haplotype listed
							  // is the bigger of the two
				firsthet=0;
				haplos[0].push_back('1');
				haplos[1].push_back('0');
				findHaplos_all(haplos,allhaps,genotype,locus+1,numloci,firsthet);
			} else {

				haplo0=haplos[0];
				haplo1=haplos[1];
				haplos[0].push_back('0');
				haplos[1].push_back('1');
				findHaplos_all(haplos,allhaps,genotype,locus+1,numloci,firsthet);

				haplos[0]=haplo0;
				haplos[1]=haplo1;
				haplos[0].push_back('1');
				haplos[1].push_back('0');
				findHaplos_all(haplos,allhaps,genotype,locus+1,numloci,firsthet);
			}
		}
	} else {
		allhaps.push_back(haplos);
	}

}

void findHaplos_approx(bool isMissing, vector<vector<string> >& allhaps,
		genoNode sample, const vector<string>& haplolist )
{
	vector<string> haplos(2);
	int total=haplolist.size();
	int numloci=sample.genotype.size();
	bool goodhaplo=0, remove=0;
	int j=0, i=0, k=0, hapindex=0;
	vector<int> indices(0);
	string genotype=sample.genotype;
	string tempstr="";

	// Make a vector of indices of haplolist
	for (i=0; i<total; i++){
		indices.push_back(i);
	}

	i=0;
	while (i<total){ // Over haplotypes in the list

		hapindex=indices[i];
		while ( (goodhaplo==0)&&(j<numloci) ){ // Over loci in haplotype if
											   // this haplo doesn't fail.
			if ( (genotype[j]=='0')&&(haplolist[hapindex][j]=='1') ){
				goodhaplo=1;
			} else if ( (genotype[j]=='2')&&(haplolist[hapindex][j]=='0') ){
				goodhaplo=1;
			}
			j++;
		}

		if (goodhaplo==0){// This haplotype can form part of a compatible pair


			haplos[0]=haplolist[hapindex];
			haplos[1]=haplos[0];

			// Find the other haplotype of the pair
			for (j=0; j<numloci; j++){
				if ( (genotype[j]!='0')&&(genotype[j]!='2') ){// Heterozygous locus
					// Find which haplo should have the 0 or 1
					if (haplolist[hapindex][j]=='0'){
						haplos[1][j]='1';
					} else {
						haplos[1][j]='0';
					}
				}
			}


			// Check if other one is in the list. If it is, remove the index.
			k=i;
			while ( (k<total)&&(remove==0) ){
				if (haplolist[indices[k]]==haplos[1]){
					remove=1;
					indices.erase(indices.begin()+k);
					total=indices.size();
				}
				k++;
			}

			// The order of the configuration matrix should be the larger haplo
			// in the first column. So swap the two if need be:

			if ((sample.hetsites.size()>0)&&(haplos[0][sample.hetsites[0]]=='0')){
				// The first one should be the larger, so swap the two
				tempstr=haplos[0];
				haplos[0]=haplos[1];
				haplos[1]=tempstr;
			}

			// Add the two haplotypes to the set of configurations
			allhaps.push_back(haplos);

		}


		// Get ready for the next round
		j=0;
		i++;
		goodhaplo=0;
		remove=0;

	}


	// TEMPORARY: Print out stuff to make sure it is doing what it should
	Rcpp::Rcout << allhaps.size() << " possible configurations for genotype: "
				 << genotype << endl;
	for (i=0; i<(int)allhaps.size(); i++){
		Rcpp::Rcout << allhaps[i][0] << " " << allhaps[i][1] << endl;
	}

	// Print out error/warning messages
	if (allhaps.size()==0) {
		/*cerr*/ Rcpp::Rcout << "ERROR: There are no haplotype configurations possible for "
				"genotype " << genotype << "." << endl;
		/*cerr*/ Rcpp::Rcout << "\tAt least one haplotype in list must be compatible with"
				" each observed genotype" << endl;
		throw Rcpp::exception("There are no haplotype configurations possible"); //exit(1);
	}

	if ( (allhaps.size()==1)&&(isMissing==1) ){
		/*cerr*/ Rcpp::Rcout << "WARNING: Only one haplotype configuration is compatible"
				" with genotype " << genotype << "." << endl;
		/*cerr*/ Rcpp::Rcout << "\t Dictionary rephaser will not be able to update this "
				" genotype" << endl;
		/*cerr*/ Rcpp::Rcout << "\t unless a compatible haplotype is added to the list of known haplotypes" << endl;
	}




}




