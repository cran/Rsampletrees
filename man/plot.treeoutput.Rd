\name{plot.treeoutput}
\alias{plot.treeoutput}
\title{Plot selected output of a sampletrees run}
\description{
This function produces the following plots:

1) Barchart giving the acceptance proportions for each update type

2) Traceplot of mutation rate (theta)

3) Traceplot of recombination rate (rho)

If tree statistics have been computed on the sampled trees, traceplots
of the tree statistics will also be produced.
}
\usage{\method{plot}{treeoutput}(x, oneperpage = FALSE, asktoplot = FALSE, layoutmat = NULL, 
	statnames = NULL, \ldots)}
\arguments{
  \item{x}{An object of class `treeoutput' containing the settings and output from a
  sampletrees run}
  \item{oneperpage}{Only one plot per page (Optional)}
  \item{asktoplot}{Plots are printed to screen after hitting enter (Optional) }
  \item{layoutmat}{A matrix giving the layout for the plots on the page. This 
  is passed to the layout() function (Optional) }
  \item{statnames}{Names for the tree statistics for plotting.}
  \item{\ldots}{Other arguments to be passed to plot.}
}
\details{
The acceptance proportions for update types 1,2, and 4-7 are the number of accepted changes 
divided by the number of times a given update was performed. Update type 3 proposes 
changes to each node of the tree in turn. The acceptance proportion for update type 3 is therefore
the average of the proportion of nodes that have an update accepted.
}
\value{
The three or more plots are plotted in the graphics window.
}
\references{Burkett KM, McNeney B, Graham J. Sampletrees and Rsampletrees: sampling gene genealogies conditional on SNP genotype data. Bioinformatics. 32:1580-2, 2016}
\author{Kelly Burkett}
\seealso{layout}
\examples{
#\dontrun{
require(ape)
mrca.age=function(tree)
{
	return(coalescent.intervals(tree)$total.depth)
}

#system.file("Examples/example_g_pars", package="Rsampletrees")
filename=paste(path.package("Rsampletrees"),"/extdata/example_g_pars",sep="")
runpars=readArgs(filename, check=FALSE)

# Must change the path so that the output can be found
#paste(system.file(package="Rsampletrees"),runpars$RunName, sep="/")
runname=paste(path.package("Rsampletrees"),"extdata",runpars$RunName,sep="/")
runpars=changeArgs(runpars, RunName=runname)

# Read in the output
results=readOutput(argob=runpars)

# Plot the output
plot(results)
plot(results, layoutmat=matrix(c(1,1,2,3),byrow=TRUE,nrow=2))
plot(results, oneperpage=TRUE, asktoplot=TRUE)

# Add a tree statistic
results=addTreeStat(results, myfunc=mrca.age, funcname="Time.to.MRCA")
plot(results)

#}
}
