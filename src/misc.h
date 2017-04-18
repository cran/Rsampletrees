
/*----------------------------------------------------------------------------

  misc.h - This is the header file for misc.cpp
  
  Kelly Burkett; Jan 2008

  --------------------------------------------------------------------------*/


#ifndef MISC_H
#define MISC_H

#include <vector>
#include <string>
#include <iostream>

using namespace std;


// This function is used to convert an integer to a string
string itos(int anumber);

// This converts an integer to a character type
char itoc(int anint);

// This converts a string to a double and string to a long long int
bool stof( const string &s, double &i);
bool stoLL( const string &s, unsigned long long int &i);


// This template function prints out a matrix formed of vector components
// The elements of the vector an be of any type. 
// NOTE: Because it is a template, the function definition is in the header
// file, though I am aware that there are work arounds. This was the only
// work around that actually worked!
template<class T>
void printMatrix(vector<vector<T> > mymatrix, ostream& outfile)
{
	int i=0, j=0;
	
	for (i=0; i<mymatrix.size(); i++){
		for (j=0; j<mymatrix[i].size()-1; j++){
			outfile << mymatrix[i][j] << " ";
		}
		outfile << mymatrix[i][j] << endl;
	}
}

// This template function prints a vector of elements of any type.  
// NOTE: Because it is a template, the function definition is in the header
// file, though I am aware that there are work arounds. This was the only
// work around that actually worked!
template<class T>
void printVector(vector<T> myvec)
{
	for (int j=0; j<myvec.size(); j++){
		cout << myvec[j] << " ";
	}
	cout << endl << endl;
}

#endif

