\name{setWeights}
\alias{setWeights}
\title{Set weights for the sampletrees proposal distributions}
\description{
Set the probabilities/weights for sampling which proposal distribution to apply in a given step
of the MCMC chain.
}
\usage{
setWeights(args, WeightFile = NULL, weightmat = NULL)
}
\arguments{
  \item{args}{
An object of class `pars' that stores arguments for
a sampletrees run
}
  \item{WeightFile}{
The name of file to write out the proposal distribution weights. 
If missing it is set to "weights"
}
  \item{weightmat}{
A matrix with the sampling weights for the proposal distributions. The first column
contains the sampling probability for the proposal distribution listed in the second column. 
If missing, default values are used. Defaults are given in the Details section.
}
}
\details{
This function is used to write a set of weights
to a file and to set the file name argument for WeightFile in the sampletrees 
settings object. The user can use default weights (if weightmat is 
NULL) or they can specify alternative weights by providing
an argument to weightmat.

The matrix weightmat must be a 2 column matrix of maximum 7 rows.
The first column contains the sampling probability for the proposal 
distribution listed in the second column. The first column must sum
to 1. The second column will consist only of some subset of the numbers
1 through 7.

The default values for the sampling weights are:

haplotype data: (0.1,0.1,0.4,0.2,0.2) for update types 1-5

genotype data: (0.1,0.1,0.35,0.1,0.1,0.125,0.125) for update types 1-7
}
\value{
\item{args}{An object of class `pars' with the tag WeightFile set to the file specified.}
	
The weights are also written out to the file specified by WeightFile.	
	}
\author{Kelly Burkett}
\references{TBD}
\examples{
runpars=newArgs()
runpars=setWeights(runpars, "weights-h.txt")

runpars=newArgs(DataType="g")
runpars=setWeights(runpars, "weights-g-2.txt", 
weightmat=cbind(c(0.05,0.05,0.4,0.2,0.15,0.15),c(1:4,6:7)))
}
