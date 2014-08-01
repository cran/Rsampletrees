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
\references{TBD}
\author{Kelly Burkett}
\examples{
settingsname=system.file("Examples/example_h_pars", package="Rsampletrees")
runpars=readArgs(settingsname, check=FALSE)
}
