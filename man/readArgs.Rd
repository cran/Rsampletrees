\name{readArgs}
\alias{readArgs}
\title{Read arguments for a sampletrees run from a settings file}
\description{
This function is used to read in settings for a sampletrees run 
that have been previously saved to a file. }
\usage{readArgs(filename, check=TRUE)}
\arguments{
  \item{filename}{The name of the input file with the sampletrees settings}
  \item{check}{If TRUE, error checking of the settings will be done (default=TRUE)}
}
\value{
Returns an object of class 'pars' with settings from the 
input file
}
\references{Burkett KM, McNeney B, Graham J. Sampletrees and Rsampletrees: sampling gene genealogies conditional on SNP genotype data. Bioinformatics. 32:1580-2, 2016}
\author{Kelly Burkett}
\examples{
filename=paste(path.package("Rsampletrees"),"/extdata/example_g_pars",sep="")
runpars=readArgs(filename, check=FALSE)
runname=paste(path.package("Rsampletrees"),"extdata",runpars$RunName, sep="/")
runpars=changeArgs(runpars, RunName=runname)
print.pars(runpars)
}
