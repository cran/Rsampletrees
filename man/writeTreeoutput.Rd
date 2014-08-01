\name{writeTreeoutput}
\alias{writeTreeoutput}
\title{
Write the results of a sampletrees run to file
}
\description{
This function can be used to write out the results of a sampletrees run to files. This function can be used when sampletrees output has been processed and added to the procdata component of the treeoutput object. It can also be used when the output of two sampletrees runs have been merged. 
}
\usage{writeTreeoutput(output, argfile=NULL, ask=TRUE)}
\arguments{
  \item{output}{An object of class 'treeoutput' with processed or merged results}
  \item{argfile}{The name for the settings file}
  \item{ask}{If the output files already exist, prompt the user to overwrite the files}
}
\details{
A settings file (argfile) and results files will be created. The name of the
results files are based on the RunName, which is stored in the runinfo component of the treeoutput object. The value of RunName is the prefix for the results files. The suffix
for the results files are '_samples.out', '_trees.out', '_accept.out' and
'_proctrees.out'.

If any of the files already exist, the user will be prompted about whether they would like to
overwrite the files or not (if ask=TRUE). If the response is `no', then the results will not be written out.
}
\references{TBD}
\author{Kelly Burkett}
\seealso{merge.treeoutput}
\examples{
mrca.age=function(tree)
{
	return(coalescent.intervals(tree)$total.depth)
}

filename=system.file("Examples/example_h_pars",package="Rsampletrees")

runpars=readArgs(filename, check=FALSE)
runname=paste(system.file(package="Rsampletrees"),runpars$RunName, sep="/")
runpars=changeArgs(runpars, RunName=runname)
results=readOutput(runpars)

results=addTreeStat(results, myfunc=mrca.age)
writeTreeoutput(results, ask=FALSE)
}