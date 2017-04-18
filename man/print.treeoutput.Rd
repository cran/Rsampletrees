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
\references{Burkett KM, McNeney B, Graham J. Sampletrees and Rsampletrees: sampling gene genealogies conditional on SNP genotype data. Bioinformatics. 32:1580-2, 2016}
\author{Kelly Burkett}
\seealso{print.pars}
\examples{
#\dontrun{
#system.file("Examples/example_g_pars",package="Rsampletrees")
filename=paste(path.package("Rsampletrees"),"/extdata/example_g_pars",sep="")

runpars=readArgs(filename, check=FALSE)
#paste(system.file(package="Rsampletrees"),runpars$RunName, sep="/")
runname=paste(path.package("Rsampletrees"),"extdata",runpars$RunName,sep="/")
runpars=changeArgs(runpars, RunName=runname)
results=readOutput(argobj=runpars)

print(results)

#}
}
