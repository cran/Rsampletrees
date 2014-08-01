\name{readTrees}
\alias{readTrees}
\title{
Read in trees from a run of sampletrees.
}
\description{
This function is used to read in trees from a sampletrees. Trees are stored as a list of
class `multiPhylo'. The user can specify that
all or a subset of the trees be read in; see Details, below. 
}
\usage{
readTrees(output=NULL, filenames=NULL,  all=TRUE, lines=NULL, 
start=1, end=NULL, nlines=NULL)
}
\arguments{
\item{output}{A treeoutput object}
\item{filenames}{A vector of names for the tree files}
  \item{all}{If TRUE, all trees in the file(s) will be read in. If FALSE, the trees specified by
  lines or start/stop/nlines will be read in}
  \item{lines}{A vector containing the line numbers of the tree file to be read in}
  \item{start}{The first line of the tree file to read in}
  \item{end}{The last line of the tree file to read in}
  \item{nlines}{The number of lines to be read in}
}
\details{
The trees are read in using the read.tree() function from the ape package, which stores results
as an object of class `multiPhylo'. 

The user has the option to provide either a vector of tree names or a treeoutput object for reading in the trees. At least one of filenames or output must not be NULL. 

To read in all of the trees in the file, use the all=TRUE option. If all=FALSE, 
then a subset of trees are read in using either the lines, start/stop or start/nlines options. 

For the lines option, the `lines' argument will consist of a vector of positive integers.  These numbers correspond to the lines in the tree file rather than to the MCMC sample number. For example, say that a chain of length 2000 is run, with a thinning interval of 100 and no burn-in. The lines of the tree file will be the 100th, 200th, 300th, etc. trees. To read in the first 4 trees, set lines=1:4 and not lines=c(100,200,300,400). This option is 
useful if non-consecutive rows are to be read in.

If the start/stop/nlines options are used, 
`start' should be set to the first row number to be read in. If a stop line is provided, 
then all lines between and including `start' and `stop' will be read in. If `nlines' is provided, 
than a total of `nlines' will be read in, starting from the row given by `start'.}
\note{Note that this function does not store the trees in a
treeoutput object; use addTrees instead.}
\value{
An object of class `multiPhylo' which is a list of tree elements.
}
\references{TBD}
\author{Kelly Burkett; uses functions from APE package, maintained by Emmanuel Paradis}
\seealso{
read.tree, addTrees
}
\examples{
treesname=system.file("Output/Example-h_trees.out", package="Rsampletrees")

# Read in all the trees in the file; may be slow
mytrees=readTrees(filename=treesname)
length(mytrees)

# Read in the first and 90th line
mytrees=readTrees(filename=treesname,all=FALSE, lines=c(1,40))
names(mytrees)

# Read in lines  2-4
mytrees=readTrees(filename=treesname,all=FALSE, lines=2:4)
names(mytrees)

# Read in three 3 starting at line 2
mytrees=readTrees(filename=treesname, all=FALSE, start=2, nlines=3)
names(mytrees)

# Read in lines 2-4
mytrees=readTrees(filename=treesname, all=FALSE, start=2, end=4)
names(mytrees)
}
