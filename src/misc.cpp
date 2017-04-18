/*----------------------------------------------------------------------------

  misc.cpp - This file contains miscellaneous functions for manipulating and
             printing data. 

  
  Notes:
  - The header file contains the functions for printing because they are
    template functions and this seemed to be the only way to reliably get
    the program to build.

  Kelly Burkett; Jan 30, 2008

  --------------------------------------------------------------------------*/


#include <string>
#include <sstream>


#include "misc.h"


// This function is used to convert an integer to a string. It first reads the
// integer into a stringstream object (stringstream provides an interface to
// manipulate strings as if they were input/output streams). The str() command
// is then used to return a copy of the string object associated with
// stringstream. For a template version of this, see
// http://www.wlug.org.nz/ConvertingAnIntegerToaStringInCpp
string itos(int anumber)
{
	stringstream ss;

	ss << anumber;
	return(ss.str());
}

// This converts an integer to a character type. It is similar to itos except
// that we must first output the contents of the stringstream to a character.
char itoc(int anint)
{
	stringstream ss;
	char mychar;

	ss << anint;
	ss >> mychar;
	return(mychar);
}




bool stof( const string &s, double &i)
{
	istringstream myStream(s);

	if (myStream>>i){
		return true;
	}
	else {
		return false;
	}
}


bool stoLL( const string &s, unsigned long long int &i)
{
	istringstream myStream(s);

	if (myStream>>i){
		return true;
	}
	else {
		return false;
	}
}






