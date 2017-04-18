\name{print.pars}
\alias{print.pars}
\title{Nice display of tags and values of a sampletrees settings object}
\description{
This function prints out the values in a sampletrees settings object in
a tabular format. The first column of the table is the names of 
all the tags. The second column gives their corresponding values.
}
\usage{\method{print}{pars}(x, \ldots)}
\arguments{
\item{x}{An object of class 'pars'}
\item{...}{additional arguments to print; currently unused}
}
\references{Burkett KM, McNeney B, Graham J. Sampletrees and Rsampletrees: sampling gene genealogies conditional on SNP genotype data. Bioinformatics. 32:1580-2, 2016}
\author{Kelly Burkett}
\seealso{summary.pars, checkArgs, print.treeoutput}
\examples{
runpars=newArgs()
print(runpars)
runpars
}

