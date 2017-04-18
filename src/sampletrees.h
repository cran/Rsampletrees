#include <iostream>
#include <ctime>
#include <fstream>
#include <cmath>
#include <sstream>
#include <bitset>
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>
#include <cstdlib> // For exit(). Newer compilers require this but it used
                   // be in iostream apparently

#include "readFiles.h"
#include "treeBuild.h"
#include "proposal.h"
#include "models.h"
#include "misc.h"
#include "rng.h"
#include "nodefunctions.h"
#include "p1p2.h"
#include "p3.h"
#include "p4.h"
#include "p5.h"
#include "p6.h"
#include "p7.h"

using namespace std; 

struct acceptOut
{
 vector<int>    one;
 vector<double> two;
 vector<double> three;
};

struct postProbOut
{
 vector<int>    i;
 vector<double> logPostProb;
};

struct RcppList
{
 acceptOut      acceptSamples;
 postProbOut    postProbSamples;
};


RcppList sampletrees ( Options user_options );
