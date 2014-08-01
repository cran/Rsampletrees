\name{print.treeoutput}
\alias{print.treeoutput}
\title{Nice display of a treeoutput object}
\description{
This function prints out information about a treeoutput object.
The component runinfo is an object of class 'pars' and is printed using print.pars. 
The total number of samples and their first/last index is also printed. }
\usage{\method{print}{treeoutput}(x, \ldots)}
\arguments{
\item{x}{An object of class 'treeoutput'}
\item{...}{additional arguments to print; currently unused}
}
\references{TBD}
\author{Kelly Burkett}
\seealso{print.pars}
\examples{
filename=system.file("Examples/example_g_pars",package="Rsampletrees")

runpars=readArgs(filename, check=FALSE)
runname=paste(system.file(package="Rsampletrees"),runpars$RunName, sep="/")
runpars=changeArgs(runpars, RunName=runname)
results=readOutput(runpars)

print(results)
}