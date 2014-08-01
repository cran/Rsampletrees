\name{addTrees}
\alias{addTrees}
\title{Read in trees and add them to the tree output object}
\description{This function will read in the trees from a from the file specified 
in the treeoutput object and store them in the treeoutput object. The user can specify that
all or a subset of the trees be read in, as described in the Details section.}
\usage{
addTrees(output, all=TRUE, lines=NULL, start=1, end=NULL, nlines=NULL)}
\arguments{
 \item{output}{An object of class `treeoutput'}
  \item{all}{If TRUE (default), all trees in the file will be read in. If FALSE, the trees specified by
  lines or start/stop/nlines will be read in}
  \item{lines}{A vector containing the line numbers of the tree file to be read in}
  \item{start}{The first line of the tree file to read in}
  \item{end}{The last line of the tree file to read in}
  \item{nlines}{The number of lines to be read in}
}
\details{ 
The trees are read in using the Rsampletrees readTree() function, which in turn calls the 
read.tree() function from the ape package. 

To read in all of the trees in the file, use the all=TRUE option. If all=FALSE, 
then a subset of trees are read in using either the lines, start/stop or start/nlines options. 

For the lines option, the `lines' argument will consist of a vector of positive integers.  These numbers correspond to the lines in the tree file rather than to the MCMC sample number. For example, say that a chain of length 2000 is run, with a thinning interval of 100 and no burn-in. The lines of the tree file will be the 100th, 200th, 300th, etc. trees. To read in the first 4 trees, set lines=1:4 and not lines=c(100,200,300,400). This option is 
useful if non-consecutive rows are to be read in.

If the start/stop/nlines options are used, 
`start' should be set to the first row number to be read in. If a stop line is provided, 
then all lines between and including `start' and `stop' will be read in. If `nlines' is provided, 
than a total of `nlines' will be read in, starting from the row given by `start'.}
\value{An object of class `treeoutput' with the trees stored in the rawdata component}
\references{TBD}
\author{Kelly Burkett}
\seealso{readOutput, readTrees}
\examples{
# Read in the settings; Must change the RunName so that the example files can be found
filename=system.file("Examples/example_g_pars",package="Rsampletrees")
runpars=readArgs(filename, check=FALSE)
runname=paste(system.file(package="Rsampletrees"),runpars$RunName, sep="/")
runpars=changeArgs(runpars, RunName=runname)


# 1. Read in all the trees; may be slow
results=readOutput(runpars)
results=addTrees(results)
length(results$rawdata$Trees)

# 2. Read in a selection of lines
results=readOutput(runpars)
results=addTrees(results, all=FALSE, lines=c(5,10,62,63,90))
length(results$rawdata$Trees)

# 3. Read in trees from lines 5 to 15
results=readOutput(runpars)
results=addTrees(results, all=FALSE, start=5, end=15)
length(results$rawdata$Trees)

# 4. Read in 20 trees, starting at line 10
results=readOutput(runpars)
results=addTrees(results, all=FALSE, start=10, nlines=20)
length(results$rawdata$Trees)
}