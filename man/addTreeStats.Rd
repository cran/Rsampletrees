\name{addTreeStat}
\alias{addTreeStat}
\title{Add tree statistics to a treeoutput object}
\description{This function computes a summary statistic on each of the trees given in the treeoutput object. The 
 results are stored in the procdata component of the treeoutput object.}
\usage{addTreeStat(output, myfunc, funcname=NULL, maxlines=1000, treerange=NULL, ...)}
\arguments{
 \item{output}{An object of class `treeoutput'}
  \item{myfunc}{The function to be applied to each tree}
  \item{funcname}{The name of the function to be used for the column name of the TreeStat data frame}
  \item{maxlines}{The maximum number of lines of the tree file to read in at a time. This avoids large tree
  files from fully being read in. }
  \item{treerange}{A vector giving the indices of the trees to which the function should be applied. Only used if
  the trees haven't been read in}
  \item{...}{further arguments that can be passed to treeapply.}
}
\details{
As with the function treeapply, the tree statistics can be computed on a subset of
all the trees (provided by treerange). Note, however, that if TreeStats already exists due to previous call to addTreeStat(), the function is applied to the 
same trees as indexed in the first column of TreeStats even if these are different from
treerange (a warning is given). This ensures a data frame of the correct dimensions.}
\value{Returns an object of class `treeoutput' that is the same as output, except that 
the tree statistics have been added to the data frame TreeStats in the procdata component.

The data frame TreeStats consists of a column for the tree index and a column for each
of the tree statistics computed. 
}
\references{Burkett KM, McNeney B, Graham J. Sampletrees and Rsampletrees: sampling gene genealogies conditional on SNP genotype data. Bioinformatics. 32:1580-2, 2016}
\author{Kelly Burkett}
\seealso{readOutput, treeapply}
\examples{
# A function that computes the time to MRCA of a tree using the ape package
require(ape)
mrca.age=function(tree)
{
	return(coalescent.intervals(tree)$total.depth)
}
#\dontrun{
 # Read in the settings; Must change the RunName so that the example files can be found
#system.file("Examples/example_g_pars",package="Rsampletrees")
 filename=paste(path.package("Rsampletrees"),"/extdata/example_g_pars",sep="")
 runpars=readArgs(filename, check=FALSE)
 
#paste(system.file(package="Rsampletrees"),runpars$RunName, sep="/")
 runname=paste(path.package("Rsampletrees"),"extdata",runpars$RunName,sep="/")
 runpars=changeArgs(runpars, RunName=runname)
 
 # Read in the output and add the summary statistic
 results=readOutput(argobj=runpars)
 results=addTreeStat(results, myfunc=mrca.age, funcname="Time.to.MRCA")
#}
}
