
#include <Rcpp.h>
#include "rcpp_sampletrees.h"
#include "sampletrees.h"
#include <string>

using namespace std;

// [[Rcpp::export]]
// inputStrings  = run_name, datatype, datafile, location_file, weightFile, HaploFreqFile, HaploListFile, InitialHaploFile, InitialTreeFile, // st_file;
// inputIntegers = seed,  ChainLength,  BurnIn,  Thinning,  RandomTree,  InitialHaplos,  HaploList,  output, InitialTree,  // st_rho, st_SA;
// inputDoubles  = FocalPoint,  InitialTheta,  MinTheta,  MaxTheta,  InitialRho,  ScaleRho,  ShapeRho;

/* This function links R with c++. The different types of data (Strings, integers and doubles) are transmitted separately by the wrappers (SEXP).
*	From SEXP, it is translated into the corresponding type vectors from the Rcpp library.
*	From the Rcpp library, it is translated into the standard library, using the Rcpp::as<> function. 
*	Example : SEXP -> Rcpp::StringVector -> std::string 
*
*  All parameters are then put into an object of class 'Options' called "user_options". This is the class that can grow to include other options such as the trios.
*	Using variables father and mother that are pointers to objects from the class 'Options'.
*  
*  The object "user_options" is transmitted to sampletrees which takes the options and returns a RcppList (see sampletrees.h) object containing :
*	acceptSamples ('acceptOut' object) that replaces 'acceptfile' variable
*	postProbSamples ('postProbOut' object) that replaces 'postprobfile' variable
*  
*  The results are sent back to R in a Rcpp::List object.
*/
RcppExport SEXP rcpp_sampletrees(SEXP inputStringsVect, SEXP inputIntegersVect, SEXP inputDoublesVect)
{
 Rcpp::StringVector  inputStrings  (inputStringsVect);
 Rcpp::IntegerVector inputIntegers (inputIntegersVect);
 Rcpp::NumericVector inputDoubles  (inputDoublesVect);
 
 Options user_options = {"",0,0,100,1,1000,"",-1,"",0,"",0,"",0,/*0,0,"",*/'h',"",-1,-1,-1,-1,-1,-1,"",0,""};

 user_options.run_name           = Rcpp::as< string >( inputStrings(0) );
 string tmp = Rcpp::as< string >( inputStrings(1) );
 user_options.datatype           = tmp.at(0);
 user_options.datafile           = Rcpp::as< string >( inputStrings[2] );
 user_options.location_file      = Rcpp::as< string >( inputStrings[3] );
 user_options.weight_file        = Rcpp::as< string >( inputStrings[4] );
 user_options.haplo_freq_file    = Rcpp::as< string >( inputStrings[5] );
 user_options.haplo_list_file    = Rcpp::as< string >( inputStrings[6] );
 user_options.initial_haplo_file = Rcpp::as< string >( inputStrings[7] );
 user_options.initial_tree_file  = Rcpp::as< string >( inputStrings[8] );
// user_options.st_file            = Rcpp::as< string >( inputStrings[9] ); // simulated tempering file
 
 user_options.seed           = inputIntegers[0];
 user_options.len_chain      = inputIntegers[1];
 user_options.burn_in        = inputIntegers[2];
 user_options.thinning       = inputIntegers[3];
 user_options.random         = inputIntegers[4];
 user_options.initial_haplos = inputIntegers[5];
 user_options.haplo_list     = inputIntegers[6];
 user_options.output         = inputIntegers[7];
 user_options.initial_tree   = inputIntegers[8];
// user_options.st_rho         = inputIntegers[9];  // simulated tempering rho 0 or 1
// user_options.st_SA          = inputIntegers[10]; // simulated tempering stochastic approximation (0 or 1)
 
 user_options.x             = inputDoubles[0]; // focalpoint
 user_options.initial_theta = inputDoubles[1];
 user_options.min_theta     = inputDoubles[2];
 user_options.max_theta     = inputDoubles[3];
 user_options.initial_rho   = inputDoubles[4];
 user_options.scale_rho     = inputDoubles[5];
 user_options.shape_rho     = inputDoubles[6];
 
 
// Rcpp::Rcout << "lancement de sampletrees" << endl;
 RcppList retour = sampletrees( user_options );
// Rcpp::Rcout << "fin de sampletrees" << endl;
 
 return Rcpp::List::create(Rcpp::Named("accept") = Rcpp::DataFrame::create( Rcpp::Named("Type")     = retour.acceptSamples.one		,
 															 Rcpp::Named("Total")    = retour.acceptSamples.two		,
 															 Rcpp::Named("Accepted") = retour.acceptSamples.three		),
					  Rcpp::Named("postProbSamples") = Rcpp::DataFrame::create(Rcpp::Named("postProb.i")    = retour.postProbSamples.i		,
 															 	    Rcpp::Named("postProb.log")  = retour.postProbSamples.logPostProb));
}

