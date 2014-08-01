\name{summary.treeoutput}
\alias{summary.treeoutput}
\title{Summary of a tree output object}
\description{Summary of a tree output object containing results of a sampletrees run}
\usage{\method{summary}{treeoutput}(object, \ldots)}
\arguments{
\item{object}{An object of class `treeoutput'}
\item{...}{additional arguments to summary; currently not used}
}
\value{Prints out a summary of the output from a sampletrees run.
The summary includes the number of MCMC samples, acceptance
proportions for each proposal type, and summary statistics (mean, median, interquartile range)
for the theta, rho and TreeStats samples.}
\references{TBD}
\author{Kelly Burkett}
\seealso{print.treeoutput, summary}
\examples{
filename=system.file("Examples/example_g_pars",package="Rsampletrees")

runpars=readArgs(filename, check=FALSE)
runname=paste(system.file(package="Rsampletrees"),runpars$RunName, sep="/")
runpars=changeArgs(runpars, RunName=runname)
results=readOutput(runpars)

summary(results)
}