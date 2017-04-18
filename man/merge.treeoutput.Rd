\name{merge.treeoutput}
\alias{merge.treeoutput}
\title{
Merge the results from two runs of sampletrees
}
\description{
This function can be used to concatenate the results of two consecutive runs of sampletrees. 
}
\usage{
\method{merge}{treeoutput}(x, y, runname=NULL, \ldots)
}

\arguments{
  \item{x}{An object of class `treeoutput' with results from the first sampletrees run}
  \item{y}{An object of class `treeoutput' with results from the second sampletrees run.}
  \item{runname}{The name to give the merged run}
  \item{\ldots}{additional arguments to merge; currently unused}
}
\details{
The function restartRun can be used to set options for starting sampletrees where 
a previous run finished. Once the second run is complete, it is useful to merge the
processed output of the two runs. This function does not merge the actual files; instead,
a new object of class `treeoutput' is created that stores the run information, the merged
Theta values, Rho values, and any tree summaries that were computed. Note that it assumes that
the same tree summaries were computed on both of the runs. To save the merged results use writeTreeoutput.

In order to set values for the first sample and last sample, the value of the 
thinning interval is used. The computation assumes that the thinning interval is the
same for both sets of output. A warning is given if they are not the same, but 
computation is not stopped; the thinning for the first set of output is assumed. 
}
\value{
Returns a list of class `treeoutput' with values set to those in output1 and output2.
}
\references{Burkett KM, McNeney B, Graham J. Sampletrees and Rsampletrees: sampling gene genealogies conditional on SNP genotype data. Bioinformatics. 32:1580-2, 2016}
\author{Kelly Burkett}
\seealso{restartRun, writeTreeoutput}
\examples{
#\dontrun{
#system.file("Examples/example_h_pars",package="Rsampletrees")
#paste(system.file(package="Rsampletrees"),runpars$RunName, sep="/")
filename=paste(path.package("Rsampletrees"),"/extdata/example_h_pars",sep="")
runpars=readArgs(filename, check=FALSE)
runname=paste(path.package("Rsampletrees"),"extdata",runpars$RunName, sep="/")
runpars=changeArgs(runpars, RunName=runname)
results1=readOutput(argobj=runpars)

#system.file("Examples/example_h_2_pars",package="Rsampletrees")
#paste(system.file(package="Rsampletrees"),runpars$RunName, sep="/")
filename=paste(path.package("Rsampletrees"),"/extdata/example_h_2_pars",sep="")
runpars=readArgs(filename, check=FALSE)
runname=paste(path.package("Rsampletrees"),"extdata",runpars$RunName, sep="/")
runpars=changeArgs(runpars, RunName=runname)
results2=readOutput(argobj=runpars)
newresults=merge(x=results1, y=results2, runname="Merge2000")

#}
}
