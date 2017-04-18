/*----------------------------------------------------------------------------

  sampletrees - This is the main file for running the MCMC algorithm to
                sample ancestral trees with internal sequence data given
                sequence data in a file.

  Kelly Burkett; March 9, 2007

  Example Usage:
  ./sampletrees optionfile runname
  where optionfile contains the user defined options
  to run the program and runname is a prefix used for all of the output files
  
  --------------------------------------------------------------------------*/


#include <gsl/gsl_rng.h>

#include "sampletrees.h"
#include <Rcpp.h>

using namespace std;


// The constants MUT_RATE and RECOMB_RATE give the per site/per pair
// mutation rate. N is used when computing theta and rho. RATE_SCALE
// values are the factor to multiply the initial value by to get the
// range of the uniform priors for theta and rho
const double RECOMB_RATE = 0.00000001;
const double MUT_RATE = 0.00000001;
const double N = 10000;
const double RATE_SCALE_MIN = 0.1;
const double RATE_SCALE_MAX = 100;


// Global variables related to the random number generator. I have
// set the random number generator to the Mercenne Twister
const gsl_rng_type * T = gsl_rng_mt19937;
gsl_rng * rng;

RcppList sampletrees ( Options user_options )
{
	unsigned int i=0, j=0;

	time_t Start_t, End_t;
	int time_task1;
	
	// Read in or set the values for the weights of the type of updates
	vector<double> weights;
	vector< vector<int> > updates;
	int dictflag=0; //Indicates that the dictionary rephaser is done

	readWeights( weights, updates, user_options.weight_file,
			     user_options.datatype, /*user_options.st_rho,*/  dictflag );


	// Set up output file prefixes, called the run name. If none is given no
	// prefix is added. Order of precedence is command line then options file.

	Rcpp::Rcout << "Output Files" << endl;

	string prefix = user_options.run_name;
	int prefix_size = prefix.size();
	
	if ( user_options.run_name!="" ){
		prefix = user_options.run_name;
		prefix +="_";
		Rcpp::Rcout << "\tOutput files start with run name: " << user_options.run_name << endl << endl;
	} else {
		Rcpp::Rcout << "\tNo run name specified" << endl << endl;
	}



	// Set up the random number generator. If a seed isn't provided, a random
	// one is used
	rng = gsl_rng_alloc (T);
	unsigned long long int myseed;

	if (user_options.seed==0){
		myseed = random_seed();
	} else {
		Rcpp::Rcout << "Using user-provided seed for RNG" << endl;
		myseed = user_options.seed;
	}
	Rcpp::Rcout << "RNG Seed: " << myseed << endl << endl;
	gsl_rng_set(rng,myseed);


	// Read in the sequence data from the sequence file specified in the
	// user options. It is read into a vector of type string.
	ifstream datafile;
	ifstream hapfile;
	ifstream inithaplofile;
	ifstream haplolistfile;
	vector<string> seq_vector;
	vector<string> haplo_list_vector(0);
	vector<genoNode> samples;
	vector<int> heterozygotes;
	vector<vector<double> > haplo_prob_mat;
	vector<vector<double> > allele_prob_mat;
	vector<vector<int> > hapindexmat;
	string temp="";
	bool found=0;

	datafile.open ( user_options.datafile.c_str() );
	if ( datafile.fail() ){
		/*cerr*/ Rcpp::Rcout << "The genotype/haplotype data file " << user_options.datafile << " cannot be opened." << endl;
		throw Rcpp::exception("The genotype/haplotype data file cannot be opened.");//return(-1); // exit ( 1 );
	}
	
	if ( user_options.datatype == 'h'){
		
		readSeq ( datafile, seq_vector );
		alleleProbs ( seq_vector, allele_prob_mat );
		haploProbs ( seq_vector, haplo_prob_mat );
		
	} else {
		
		// Read in the list of likely haplotypes. These will be passed to
		// readGenos
		if (user_options.haplo_list==1){
			haplolistfile.open ( user_options.haplo_list_file.c_str() );
			if ( haplolistfile.fail() ){
				/*cerr*/ Rcpp::Rcout << "The file listing likely haplotypes, " << user_options.haplo_list_file << ", cannot be opened." << endl;
				throw Rcpp::exception("The file listing likely haplotypes, cannot be opened."); //return(-1); // exit ( 1 );
			} else {
				readSeq( haplolistfile, haplo_list_vector);
			}
			haplolistfile.close();

		}


		// Read in the genotypes from the datafile and set up the vector with sample information
		readGenos( datafile, samples, heterozygotes, hapindexmat, dictflag, haplo_list_vector);
		
		
		// Read in the initial haplotype frequencies since this currently is computed by the program
		hapfile.open( user_options.haplo_freq_file.c_str() );
		if ( hapfile.fail() ){
			/*cerr*/ Rcpp::Rcout << "The haplotype frequency file  "
				 << user_options.haplo_freq_file
				 << "cannot be opened." << endl;
			throw Rcpp::exception("The haplotype frequency file cannot be opened."); //return(-1); // exit ( 1 );
		}
		readHaps( hapfile, haplo_prob_mat, (samples[0].genotype).size()-1 );
		hapfile.close();


		// Determine an initial haplotype configuration.
		if ( user_options.initial_haplos == 0 ){
			// No intial haplotypes were provided
			if ( (dictflag==1)&&(user_options.haplo_list==1) ){
				// Dictionary rephaser will be used and a haplotype list was
				// provided. Take initial haplotypes from enumeration
				initialHaplos_enum( seq_vector, samples );
			} else {
				//use the haplotype frequencies to determine an initial
				// phase configuration.
				initialHaplos_freqs( seq_vector, samples, haplo_prob_mat );
			}
		} else {
			inithaplofile.open ( user_options.initial_haplo_file.c_str() );
			if ( inithaplofile.fail() ){
				/*cerr*/ Rcpp::Rcout << "The initial haplotype data file "
					 << user_options.initial_haplo_file
				     << "cannot be opened." << endl;
				throw Rcpp::exception("The initial haplotype data file cannot be opened."); //return(-1); // exit ( 1 );
			}
			readSeq ( inithaplofile, seq_vector );

			// Set up the initial sequence ids
			for (i=0; i<samples.size(); i++){
				if ((samples[i].hetsites.size()!=0)&&(seq_vector[2*i][samples[i].hetsites[0]]=='0')){
					// Swap the two sequences so that the largest in decimal
					// numbers is the first sequence
					temp=seq_vector[2*i];
					seq_vector[2*i]=seq_vector[2*i+1];
					seq_vector[2*i+1]=temp;
				}
				samples[i].seq1 = i*2;
				samples[i].seq2 = i*2+1;
			}
			i=0;

			// Error check that if the dictionary is used and a haplo list was
			// provided that the initial haplotypes provided are actually
			// possible given the haplotype list
			if ( (dictflag==1)&&(user_options.haplo_list==1) ){
				for (i=0; i<samples.size(); i++){
					j=0;
					found=0;
					while ( (j<samples[i].hapconfigs.size())&&(found==0) ){
						if ( (samples[i].hapconfigs[j][0]==seq_vector[samples[i].seq1])||
							 (samples[i].hapconfigs[j][0]==seq_vector[samples[i].seq2]) ){
							found=1;
						}
						j++;
					}
					if (found==0){
						/*cerr*/ Rcpp::Rcout << "ERROR: Initial haplotypes provided for individual "
							 << i+1 << " are not in haplotype list" << endl;
						throw Rcpp::exception("ERROR: Initial haplotypes provided for individual are not in haplotype list"); //exit(1);
					}
				}
			}

		}

		// You can use the seq_vector to find allele_prob_mat because the number
		// of alleles is the same regardless of the sequence configuration:
		alleleProbs ( seq_vector, allele_prob_mat );
		
	}
	printSeq ( seq_vector );
	datafile.close();


	// Set up for the dictionary rephaser
	vector<int> currenthapvec(samples.size(),0);
	vector<double> P9scores(hapindexmat.size(),0);


	// Read in the marker location data from the file specified in the user
	// options. It is read into a vector
	ifstream distfile;
	vector<double> dist_vector ( seq_vector[0].length() );
	double diff=0;

	distfile.open ( user_options.location_file.c_str() );
	if ( distfile.fail() ){
		/*cerr*/ Rcpp::Rcout << "The marker location file " << user_options.location_file;
		/*cerr*/ Rcpp::Rcout << " cannot be opened." << endl;
		throw Rcpp::exception("The marker location file cannot be opened."); //return(-1); // exit ( 1 );
	}

	readDist ( distfile, dist_vector );
	distfile.close();
	diff = dist_vector[dist_vector.size()-1]-dist_vector[0];


	// Set up the initial value and the range of theta. These values may be
	// given by the user but not necessarily.
	double theta=4*N*MUT_RATE*diff;
	double theta_min=theta*RATE_SCALE_MIN;
	double theta_max=theta*RATE_SCALE_MAX;

	if (user_options.min_theta!=-1){ theta_min = user_options.min_theta; }

	if (user_options.max_theta!=-1){ theta_max = user_options.max_theta; }

	if (user_options.initial_theta!=-1){
		theta = user_options.initial_theta;
	}
	else {
		// Uniform range was provided so choose the midpoint as the initial
		// theta
		if ( (user_options.min_theta!=-1)||(user_options.max_theta!=-1) ){
			theta=(theta_min+theta_max)/2;
		}
	}

	if ( (theta<theta_min)||(theta>theta_max) ) {
		/*cerr*/ Rcpp::Rcout << "Initial Theta is outside of the bounds of the prior" << endl;
		throw Rcpp::exception("Initial Theta is outside of the bounds of the prior"); //exit(1);
	}
	Rcpp::Rcout << "Theta/Rho Values" << endl;
	Rcpp::Rcout << "\t\tInitial Theta: " << theta << endl;
	Rcpp::Rcout << "\t\tMinimum Theta: " << theta_min << endl;
	Rcpp::Rcout << "\t\tMaximum Theta: " << theta_max << endl << endl;


	// Set up the initial value and the range of rho. These values may be
	// given by the user but not necessarily.
	double rho=0;
	double rho_shape=1;
	double rho_scale=0.1;

	vector<double> normc(0);
	vector<double> lambdas(0);
	ifstream stfile;
	unsigned int I=1;

	string sttype="shape"; //shape or scale

	if (user_options.shape_rho!=-1){
		rho_shape = user_options.shape_rho;
	}

	if (user_options.scale_rho!=-1){
		rho_scale = user_options.scale_rho;
	}
/*
	if (user_options.st_rho==1){

		// Going to do simulated tempering. In this case, read in the
		// distributional information from the file specified
		stfile.open ( user_options.st_file.c_str() );
		if ( stfile.fail() ){
			 Rcpp::Rcout << "The simulated tempering file " << user_options.st_file
				 << " cannot be opened." << endl;
			throw Rcpp::exception("The simulated tempering file cannot be opened."); //return(-1); // exit ( 1 );
		}

		readst(normc, lambdas, stfile, user_options);
		I=normc.size();

		rho=lambdas[lambdas.size()-1];

		Rcpp::Rcout << "\t\tInitial Rho: " << rho << endl;
		Rcpp::Rcout << "\t\tTemperature: " << sttype << endl;
		Rcpp::Rcout << "\t\tSimulated Tempering pseudoprior and lambda values:" << endl;

		for (i=0; i<normc.size(); i++){
			Rcpp::Rcout << "\t\t\t" << normc[i] << " " << lambdas[i] << endl;
		}
		Rcpp::Rcout << endl;

	} else {

		if (user_options.initial_rho!=-1){
			rho = user_options.initial_rho;}
		else {
			rho = rho_shape*rho_scale;
		}

		Rcpp::Rcout << "\t\tInitial Rho: " << rho << endl;
		Rcpp::Rcout << "\t\tRho Shape parameter: " << rho_shape << endl;
		Rcpp::Rcout << "\t\tRho Scale parameter: " << rho_scale << endl << endl;

	}
*/



	// Compute the marker indices of the first markers to the left and right
	// of the focal point. LR_index[0] is the first marker to the left of the
	// focal point and LR_index[1] is the first marker to the right
	vector<int> LR_index ( 2,0 );
	i=0;
	
	while ( (dist_vector[i] < user_options.x) && (i < dist_vector.size()) ){
	// Find the marker that is immediately before the reference point
		i++;
	}

	
	LR_index[0]=i-1;
	LR_index[1]=i;

	if ( dist_vector[i] == user_options.x ){ 
		//the focal pt is one of the markers
		LR_index[1]=i+1;
	}


	// Set up the initial tree structure and data. We can do this either from a
	// random initial point or from using an initial data structure that was
	// specified in the options file. This will typically be output from a
	// previous run and will be of use for continuing a previously started chain.
	treeNode *head_node;

	if ( user_options.initial_tree == 0 ){
		// No initial tree data was given in the options file so generate it
		// If a random tree is desired, the sequences are not coalesced first
		// based on similarity, they are just done randomly. Otherwise make
		// a sensible tree using UPGMA
		if (user_options.random == 1){
			head_node = createRandomTree ( seq_vector, LR_index,
					                       dist_vector, user_options.x, rho);
		} else {
			head_node = createInitialTree ( seq_vector, LR_index,
											dist_vector, user_options.x, rho);
		}
	} else {
		// The initial tree data was given in the options file so use it
		// for starting conditions.
		ifstream initial_file;
		initial_file.open ( user_options.initial_tree_file.c_str() );

		if ( initial_file.fail() ){
			/*cerr*/ Rcpp::Rcout << "The initial tree data file ";
			/*cerr*/ Rcpp::Rcout << user_options.initial_tree_file << " cannot be opened.";
			/*cerr*/ Rcpp::Rcout << endl;
			throw Rcpp::exception("The initial tree data file cannot be opened."); //return(-1); // exit ( 1 );
		}
		head_node = readInitialTree ( initial_file );
		initial_file.close();
	}

	
	// Print out the full initial data to a file
	ofstream initial_tree;
	initial_tree.open( prefix.append("firsttree.out").c_str() );  // The first tree is printed out in a file if the user wants it.
	printTreeFull ( head_node, initial_tree );
	initial_tree.close();


	// Get a vector of pointers to each of the tree nodes. This will be used
	// for sampling nodes from the tree and is generally useful to have to
	// access nodes without traversing the tree. If we have genotype data,
	// relate the genotype data vector to the tree nodes. 
	vector<treeNode*> node_ptrs ( 2*seq_vector.size()-1 );
	getNodePtrs ( node_ptrs, head_node );


	// Sample from the 6 proposal distributions according to the probabilities
	// specified either in the options file or the default values.
	int proposal_index=0, accept_flag=0, numintnodes=seq_vector.size()-1;
	vector<double> temp_vec(8,0);
	vector<vector<double> > accept_mat(3,temp_vec);
//	int m=normc.size();
//	double c=1, n=123;
	vector<string> newseqs(2);

	for (i=0; i<8; i++){ accept_mat[0][i]=i+1; }
	
	ofstream statfile;
	statfile.open ( prefix.replace((prefix_size+1),(prefix.size()-prefix_size-1),"samples.out").c_str() );
	statfile << "i theta rho" << endl;
	
	ofstream treefile;
	treefile.open ( prefix.replace((prefix_size+1),(prefix.size()-prefix_size-1),"trees.out").c_str() );
	
//	ofstream acceptfile;
//	string acceptfilename= prefix.replace((prefix_size+1),(prefix.size()-prefix_size-1),"accept.out").c_str();
	acceptOut acceptSample ;		// <- replacement
	
	//ofstream postprobfile;
	postProbOut postProbSample;	// <- replacement
/*
	if (user_options.st_rho==1){
		stoutfile.open ( prefix.replace(prefix_size,11,"st.out").c_str() );
		Stout.push_back("I");
		stoutfile << "I"<< endl;
		rfile.open ( prefix.replace(prefix_size,10,"rs.out").c_str() );
	}
*/
	ofstream samplefile;
	
	ofstream haplooutfile;
	string haplolistname;

	if (user_options.datatype=='g'){
		samplefile.open( prefix.replace((prefix_size+1),(prefix.size()-prefix_size-1),"sequences.out").c_str() );
		if (user_options.haplo_list==1){
			haplolistname=prefix.replace((prefix_size+1),(prefix.size()-prefix_size-1),"haplolist.out");
		}
	}


	// Return the log of the posterior probability as a convergence diagnostic
	//postprobfile.open( prefix.replace(prefix_size,13,"postprob.out").c_str() );


	// Set up for the output file for the last tree. We will overwrite this file just
	// in case the job is stopped for whatever reason
	ofstream last_tree;
	prefix.replace((prefix_size+1),(prefix.size()-prefix_size-1),"lasttree.out");
	int update_type=9, num_updates=5;
	bool updateflag=1;

	Start_t = time(NULL);

/*	if ( user_options.st_rho==1 ){
			num_updates+=3;
	} else 
*/	if (user_options.datatype=='g'){
		num_updates++;
		if (dictflag==1){
			num_updates++;
		}
	}

//	Rcpp::Rcout << "Début de boucle "<<endl;
	for ( i=1; i<user_options.len_chain+1; i++ )
	{
		proposal_index = samplePMF ( weights );

		for ( j=0; j< updates[proposal_index].size(); j++)
		{
//			Rcpp::Rcout << "boucle " << i << " - " << j << endl;

			update_type = updates[proposal_index][j];
//			Rcpp::Rcout << "update " << update_type  << endl;
			if ( update_type == 1 ){
				accept_flag = updateP1 ( theta, node_ptrs, LR_index, theta_min, theta_max );
				accept_mat[1][0]++;
				accept_mat[2][0]+=accept_flag;
			}
			else if ( update_type == 2 ){
//				if (user_options.st_rho==0){
				accept_flag = updateP2_gamma ( rho, node_ptrs, LR_index, dist_vector,
											   user_options.x, rho_shape, rho_scale );
				accept_mat[1][1]++;
				accept_mat[2][1]+=accept_flag;
//				Rcpp::Rcout << "update 2 accept " <<  accept_flag << endl;
//				}
			}
			else if ( update_type == 3 ){
				accept_flag = updateP3 ( node_ptrs, LR_index,user_options.x,
										 dist_vector, rho, theta, allele_prob_mat,
			                             haplo_prob_mat, i );
				accept_mat[1][2]++;

				// Store the average of the proportion accepted
				accept_mat[2][2]=((accept_mat[1][2]-1)*accept_mat[2][2]+
						(double(accept_flag)/numintnodes))/accept_mat[1][2];
			}
			else if ( update_type == 4 ){
					accept_flag = updateP4 ( node_ptrs, theta, rho, LR_index,
							user_options.x,dist_vector,haplo_prob_mat,
							allele_prob_mat );
					accept_mat[1][3]++;
					accept_mat[2][3]+=accept_flag;
			}
			else if ( update_type == 5 ){
				accept_flag = updateP5 ( node_ptrs, theta, rho, LR_index,
			                             user_options.x, dist_vector, haplo_prob_mat,
			                             allele_prob_mat );
				accept_mat[1][4]++;
				accept_mat[2][4]+=accept_flag;
			}
			else if ( update_type == 6 ) {

				if (user_options.datatype!='g'){
					Rcpp::Rcout << "Can't perform this update as it is for genotype data" << endl;
				} else {
					accept_flag = updateP7_sequence( node_ptrs, samples, heterozygotes,
							haplo_prob_mat, allele_prob_mat,
							theta, LR_index, "swap", 1, updateflag, newseqs );

					if ( (accept_flag==1)&&(user_options.haplo_list==1) ){
						updateHaploList( newseqs, haplo_list_vector, samples,
								hapindexmat, P9scores );
					}
					accept_mat[1][5]++;
					accept_mat[2][5]+=accept_flag;
				}
			}
			else if (update_type == 7){
				if (user_options.datatype!='g'){
					Rcpp::Rcout << "Can't perform this update as it is for genotype data" << endl;
				} else {
					accept_flag = updateP7_topo ( node_ptrs, samples, heterozygotes,
							haplo_prob_mat, allele_prob_mat,
							theta, rho, LR_index, dist_vector, user_options.x,
							P9scores, hapindexmat, currenthapvec, updateflag );
					accept_mat[1][6]++;
					accept_mat[2][6]+=accept_flag;
				}
			}
			else if ( update_type == 8 ){
				//stringstream ssStout;
				updateI_rho( I, normc, lambdas, rho, /*stoutfile, ssStout,*/ node_ptrs, LR_index, dist_vector,	user_options.x);
				//Stout.push_back( ssStout.str() );
				// Use stochastic approximation to change the weights if desired
/*				if (user_options.st_SA==1){
					for (j=1; j<=normc.size(); j++){
						if (j==I){
							normc[j-1] = exp( log(normc[j-1])-c/(i+n) );}
						else {
							normc[j-1] = exp( log(normc[j-1])+c/(m*(i+n)) );}
					}
				} */
				accept_mat[1][7]++;

			}
			else {
				/*cerr*/ Rcpp::Rcout << "ERROR: An invalid proposal type was picked." << endl;
				throw Rcpp::exception("ERROR: An invalid proposal type was picked."); //return(-1); // exit ( 1 );
			}

		}

		if ( (i >= user_options.burn_in) && ( (i % user_options.thinning) == 0) ){

			statfile << i << " " << theta << " " << rho << endl;

			printTreeStruc ( node_ptrs[node_ptrs.size()-1], i, 1, treefile );

			if (user_options.datatype=='g'){
				samplefile << i << " ";
				//samSample.i.push_back(i);
				//stringstream samTmp;
				for (j=0; j<seq_vector.size(); j++){
					samplefile << node_ptrs[j]->seq << " ";
					//samTmp << node_ptrs[j]->seq << " ";
				}
				samplefile << endl;
				//samSample.seq.push_back( samTmp.str() );
			}

			/*postprobfile <<  i << " " << probAll_relative(node_ptrs, rho,
					theta, LR_index, dist_vector, user_options.x, allele_prob_mat,
					haplo_prob_mat, user_options.st_rho, rho_shape, rho_scale) << endl;*/
			postProbSample.i.push_back(i);
			
			postProbSample.logPostProb.push_back( 
				probAll_relative(node_ptrs, rho, theta, LR_index, dist_vector, 
					user_options.x, allele_prob_mat, haplo_prob_mat, 
					rho_shape, rho_scale));

//			acceptfile.open ( acceptfilename.c_str(), ios::trunc );
			acceptSample = {}; 
			for (int l=0; l<num_updates; l++){
				acceptSample.one.push_back  (accept_mat[0][l]);
				acceptSample.two.push_back  (accept_mat[1][l]);
				acceptSample.three.push_back(accept_mat[2][l]);
/*				for (j=0; j<3; j++){
					acceptfile << accept_mat[j][l] << " ";
				}
				acceptfile << endl;
*/			}
//			acceptfile.close();
		}

		if ( ( i % 1000 ) == 0 ){

			// Save a backup version of the tree to start the run again
			// if it is stopped.
			last_tree.open ( prefix.c_str(), ios::trunc );
			printTreeFull ( node_ptrs[node_ptrs.size()-1], last_tree);
			last_tree.close();


			// Print out the final haplotype list if it was an option
			haplooutfile.open( haplolistname.c_str(), ios::trunc );
			if ( (user_options.datatype=='g')&&(user_options.haplo_list==1) ){
				for (unsigned int l=0; l<haplo_list_vector.size(); l++){
					haplooutfile << haplo_list_vector[l] << endl;
				}
			}
			haplooutfile.close();

//			acceptfile.open ( acceptfilename.c_str(), ios::trunc );
			acceptSample = {}; 
			for (int l=0; l<num_updates; l++){
				acceptSample.one.push_back  (accept_mat[0][l]);
				acceptSample.two.push_back  (accept_mat[1][l]);
				acceptSample.three.push_back(accept_mat[2][l]);
/*				for (j=0; j<3; j++){
					acceptfile << accept_mat[j][l] << " ";
				}
				acceptfile << endl;
*/			}
//			acceptfile.close();


			// Print out information about current iteration
			if ( user_options.output == 0 ){
				Rcpp::Rcout << "Iteration: " << i << endl;
			} else {
				Rcpp::Rcout << "Iteration: " << i << " " << endl;
				printTreeFull ( node_ptrs[node_ptrs.size()-1], Rcpp::Rcout );
				Rcpp::Rcout << endl;
			}
		}
	}


	// MCMC cycles are complete.

	// Return information about the length of the run
	End_t = time(NULL);    //record time that task 1 ends
	time_task1 = difftime(End_t, Start_t);

	Rcpp::Rcout << "\n\nCompleted iterations in " << time_task1 << " seconds" << endl;


	// Close the output files
	statfile.close();
	treefile.close();
//	acceptfile.close();
	//postprobfile.close();
	if (user_options.datatype=='g'){
		samplefile.close();
	}


	// Print out the final haplotype list if it was an option
	haplooutfile.open( haplolistname.c_str(), ios::trunc );
	if ( (user_options.datatype=='g')&&(user_options.haplo_list==1) ){
		for (i=0; i<haplo_list_vector.size(); i++){
			haplooutfile << haplo_list_vector[i] << endl;
		}
	}
	haplooutfile.close();

	// Print out the full data for the final tree sampled
	last_tree.open ( prefix.c_str(), ios::trunc );
	printTreeFull ( node_ptrs[node_ptrs.size()-1], last_tree );
	last_tree.close();
	
	// Clean up
	deleteTree ( node_ptrs[node_ptrs.size()-1] );
	gsl_rng_free(rng);

	RcppList returnList;
	returnList.acceptSamples = acceptSample;
	returnList.postProbSamples = postProbSample;
	
	return returnList;
}

