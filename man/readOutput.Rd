\name{readOutput}
\alias{readOutput}
\title{Read in sampletrees output}
\description{This function is used to read in results from a run of sampletrees. 
Either the object with the settings for the run  or the
file containing these settings is passed to the function. A new
object is created that contains a component that stores the run settings, plus additional
components storing the output, the sampled values for theta and rho, and a matrix with acceptance
proportions. Note that the trees themselves are by default not read in as these files can 
be very large.  
}
\usage{
readOutput(argobj=NULL, argfile=NULL, addtrees=FALSE)
}
\arguments{
\item{argobj}{A settings object of class `pars' describing the sampletrees run}
\item{argfile}{A file containing the settings describing the sampletrees run}
\item{addtrees}{If TRUE, store the trees in the output object (default=FALSE)}
}
\details{One of either argobj or argfile must not be NULL.}
\value{
An object of class `treeoutput', which is a list made up of three components:

1) runinfo - a copy of argobj

2) rawdata

3) procdata

The component `rawdata' consists of
\item{i}{Iteration numbers of the MCMC samples}
\item{Theta}{A vector of sampled theta values (mutation rate)}
\item{Rho}{A vector of sampled rho values (recombination rate)}
\item{Trees}{Either a string containing the name of the tree file or a list
of class `multiPhylo' containing the trees (if addtrees=TRUE). See the ape package
documentation for more information on the `multiPhylo' class.}

The component `procdata' (processed data) is also a list, initially made up only of
a matrix with the acceptance proportions for each update type. Tree statistics may be added to procdata by 
addTreeStat}
\references{TBD}
\author{Kelly Burkett}
\seealso{addTreeStat, addTrees}
\examples{
filename=system.file("Examples/example_g_pars", package="Rsampletrees")
runpars=readArgs(filename, check=FALSE)

# Must change the path so that the output can be found
runname=paste(system.file(package="Rsampletrees"),runpars$RunName, sep="/")
runpars=changeArgs(runpars, RunName=runname)

# Read in the output
results=readOutput(argobj=runpars)
}
