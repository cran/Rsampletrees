\name{checkArgs}
\alias{checkArgs}
\title{
Error checking for the arguments to sampletrees
}
\description{
This function is used to check the arguments to be used for a run
of sampletrees. If the options are properly specified, then sampletrees
can be expected to run without errors from improper input files. 
}
\usage{checkArgs(args)}
\arguments{
  \item{args}{An object of class `pars' with the arguments for the sampletrees run}
}
\details{
Extensive checking is done, including:

- Non-default values have been provided for settings requiring user-specified values

- All input files exist

- All numeric settings are set to numeric values

- The format of the genotype and/or sequence files are correct

- The file with SNP locations contains only numeric values in increasing order, and 
that the number of locations matches the number of SNPs in the genotype/sequence file

- Focal point is set to a location in the region containing the SNPs

- Weights sum to 1 and the proposal indices are in the right range
}
\value{
Returns an object of class 'pars' with the 'clean' variable set to TRUE
if the arguments in `args' pass the error check and FALSE otherwise.
}
\references{
Burkett KM, McNeney B, Graham J. Sampletrees and Rsampletrees: sampling gene genealogies conditional on SNP genotype data. Bioinformatics. 32:1580-2, 2016
}
\author{Kelly Burkett}
\examples{

# Will produce error messages
runpars=newArgs()
runpars=checkArgs(runpars)

# Will not produce errors
#system.file("Examples/sequences_Theta8_Rho8.txt",package="Rsampletrees")
datname=paste(path.package("Rsampletrees"),"/extdata/sequences_Theta8_Rho8.txt",sep="")
#system.file("Examples/locations_Theta8_Rho8.txt",package="Rsampletrees")
locname=paste(path.package("Rsampletrees"),"/extdata/locations_Theta8_Rho8.txt",sep="")
#system.file("Examples/weights-h.txt", package="Rsampletrees")
weightname=paste(path.package("Rsampletrees"),"/extdata/weights-h.txt",sep="")
runpars=newArgs(DataFile=datname, DataType="h", WeightFile=weightname,
		LocationFile=locname,RunName="Test-h",FocalPoint=10000)
runpars=checkArgs(runpars)

}

