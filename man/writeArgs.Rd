\name{writeArgs}
\alias{writeArgs}
\title{
Write arguments for a sampletrees run to a file in the format required
for sampletrees
}
\description{
This function is used to write the arguments for a sampletrees run to a file. 
The arguments are stored in the settings object. The settings are written out
in the format required by sampletrees. The name of the settings file can then be used as an
argument for the run of sampletrees. 
}
\usage{
writeArgs(args, outfile = NULL)
}
\arguments{
  \item{args}{A settings object of class `pars' containing settings for sampletrees}
  \item{outfile}{The name of the output file for writing the settings}
}
\details{A warning is given if the settings to be written either 
haven't been checked using checkArgs() or if there is an error. This
warning will be given if the list value `clean' is FALSE. However, 
the settings will be written to the file indicated even with this warning.}
\references{TBD}
\author{Kelly Burkett}
\seealso{checkArgs()}
\examples{
## Note that this example will give a warning message about the settings 
## not being checked. 

runpars=newArgs(DataFile="sequences_Theta8_Rho8.txt", 
LocationFile="locations_Theta8_Rho8.txt",
RunName="Test-h",FocalPoint=10000, WeightFile="weights-h.txt")
writeArgs(runpars, "test_h_pars")
}
