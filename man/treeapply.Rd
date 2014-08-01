\name{treeapply}
\alias{treeapply}
\title{Apply a function to the sampled trees}
\description{This function is used to apply a user-supplied function to each tree that
was sampled in a run of sampletrees.  }
\usage{
treeapply(output, myfunc, funcname=NULL, maxlines=1000, treerange=NULL)
}
\arguments{
  \item{output}{An object of class `treeoutput'}
  \item{myfunc}{The function to be applied to each tree}
  \item{funcname}{The name of the function to be used for the column name of the data frame that is returned}
  \item{maxlines}{The maximum number of lines of the tree file to read in at a time. This avoids large tree
  files from fully being read in. }
  \item{treerange}{A vector giving the indices of the trees to which the function should be applied. Only used if
  the trees are not stored in the treeoutput object}
}
\details{This function is meant to be similar to the various `apply' functions, such as
lapply and sapply, in that it can be used to apply a function to each elements of a list
or vector. In this case, it applies the specified function to an object of class 'multiPhylo', 
which is a list of trees. 

The trees to apply the function to are found in the rawdata component of the treeoutput object (either the tree file name or the trees stored as an object of class `multiPhylo'). 

The function 
must be written to manipulate tree data of type 'phylo'. See the 'ape' package documentation 
for more information about the format of tree data. }
\value{Returns a data frame consisting of a column for the tree index and a column with the values of the function computed on the corresponding tree}
\references{TBD}
\author{Kelly Burkett}
\seealso{addTreeStats}
\examples{
mrca.age=function(tree)
{
	return(coalescent.intervals(tree)$total.depth)
}

filename=system.file("Examples/example_g_pars", package="Rsampletrees")
runpars=readArgs(filename, check=FALSE)

# Must change the path so that the output can be found
runname=paste(system.file(package="Rsampletrees"),runpars$RunName, sep="/")
runpars=changeArgs(runpars, RunName=runname)

# Read in the output
results=readOutput(argob=runpars, addtrees=TRUE)
test=treeapply(results, myfunc=mrca.age,funcname="TimeMRCA")

results=readOutput(argob=runpars)
test=treeapply(results, myfunc=mrca.age, treerange=50:60)

}
