\name{launch.sampletrees}
\alias{launch.sampletrees}
\docType{package}
\title{
Launches the c++ sampletrees program using arguments previously given and checked.
}
\description{
This function performs the MCMC tree sampling by launching the c++ sampletrees program. 
Linking of the pre_sampletrees functions to the c++ sampletrees program is done using Rcpp package.
For more information about the sampletrees program, please see the package vignette.
}
\usage{launch.sampletrees(args, addtrees=FALSE)}
\arguments{
  \item{args}{An object of class `pars' with the arguments for the sampletrees run}
  \item{addtrees}{If TRUE, store the trees in the output object (default=FALSE)}
}
\value{
An object of class `treeoutput', which is a list made up of three components:

1) runinfo - a copy of argobj

2) rawdata

3) procdata

The component `rawdata' consists of
\item{i}{Iteration numbers of the MCMC samples}
\item{Theta}{A vector of sampled theta values (mutation rate)}
\item{Rho}{A vector of sampled rho values (recombination rate)}
\item{Trees}{Either a string containing the name of the tree file or a list
of class `multiPhylo' containing the trees (if addtrees=TRUE). See the ape package
documentation for more information on the `multiPhylo' class.}

The component `procdata' (processed data) is also a list, initially made up only of
a matrix with the acceptance proportions for each update type. Tree statistics may be added to procdata by 
addTreeStat}
\details{
Sampletrees is a computationally intensive program. As such, it will write important MCMC output 
to output files as a back-up. Upon completion of the tree sampling, launch.sampletrees will read
the output into R as a treeoutput object. 
}
\references{Burkett KM, McNeney B, Graham J. Sampletrees and Rsampletrees: sampling gene genealogies conditional on SNP genotype data. Bioinformatics. 32:1580-2, 2016}
\examples{
\dontrun{
oldPath=getwd()
setwd(paste(path.package("Rsampletrees"),"/extdata/",sep=""))
runpars=readArgs("example_h_pars")
runtree=launch.sampletrees(runpars)
setwd(oldPath)
}
}
