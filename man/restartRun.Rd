\name{restartRun}
\alias{restartRun}
\title{
Set up a settings object for continuing a previous sampletrees run
}
\description{
This function is used to initialize settings in order to restart a 
sampletrees run. The initial values for theta, rho and the initial tree are
taken from the final sampled values of a previous sampletrees run.  }
\usage{
restartRun(newrunname, oldargs = NULL, argfile = NULL, extrait = NULL, 
totalsamples = NULL)
}
\arguments{
  \item{newrunname}{The name to associate with the new run. Output files will have
  this run name as prefix}
  \item{oldargs}{An object of class `pars' with the settings for the previous
  sampletrees run}
  \item{argfile}{The name of the settings file used for the
  previous sampletrees run}
  \item{extrait}{The number of additional iterations desired}
  \item{totalsamples}{The total number of iterations desired in the previous and new run}
}
\details{
The settings of the previous sampletrees runs can be specified in terms
of the settings object or a file name. Therefore, at least one of `oldargs' 
or 'argfile' must not be NULL.

The desired number of MCMC samples can be specified either in terms of the
additional iterations to run ('extrait') or in terms of the total number of 
iterations desired in both the previous and new run ('totalsamples'). Therefore, 
at least one of 'extrait'  or 'totalsamples' must not be NULL.

The settings in the new settings object are the same as the previous except:

1) The initial theta value is set to the last sampled value from the previous run

2) The initial rho value is set to the last sampled value from the previous run

3) The data for the initial tree, including the node times, internal sequence and
recombination-related latent variables, are set to the last sampled values from
the previous run. 

}
\value{
Returns an object of class `pars' with the settings for a sampletrees run 
that starts where the previous run finished.
}
\references{TBD}
\author{Kelly Burkett}
\examples{

filename=system.file("Examples/example_h_pars",package="Rsampletrees")
runpars=readArgs(filename, check=FALSE)

# Include path in run name so that function can find the necessary files
runname=paste(system.file(package="Rsampletrees"),runpars$RunName, sep="/")
runpars=changeArgs(runpars,RunName=runname)

newpars=restartRun("example-h-2.pars", oldargs=runpars, totalsamples=200000)
}