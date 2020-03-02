\name{estimateHap}
\alias{estimateHap}
\title{Estimate and write haplotype probabilities and initial values to files}
\description{
This function is only used when the data type is `g' (genotype).
Estimate the two-locus haplotype probabilities from the genotype data. 
Optionally set initial haplotype configurations and a list of likely 
haplotypes for the run of sampletrees. This function
requires the R package "haplo.stats" for estimating the haplotype 
frequencies from the genotype data.
}
\usage{
estimateHap(args, HaploFreqFile, InitialHaplos = TRUE, 
InitialHaploFile = "initialhaps", HaploList = TRUE, 
HaploListFile = "initialhaplist", tol = 1e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{args}{An object of class `pars' with the sampletrees settings}
  \item{HaploFreqFile}{The name of the output file for the two-locus haplotype 
  probability estimates}
  \item{InitialHaplos}{If TRUE, sample the initial haplotype configuration
  using the posterior probabilities of each configuration for each individual
  estimated using haplo.em() (Default=TRUE)}
  \item{InitialHaploFile}{File name for the initial haplotype
  configurations}
  \item{HaploList}{If TRUE, create a list of likely haplotypes for a run
  of sampletrees with genotype data (Default=TRUE). This list will include haplotypes estimated
  to have a probability greater than 'tol'}
  \item{HaploListFile}{File name for the haplotype list}
  \item{tol}{Haplotypes with estimated probability greater than this value will be included in the list of likely 
  haplotypes}
}
\details{
This function is only used when the data type is genotype (`g').

The two-locus haplotype probabilities are estimated using the haplo.em() function in the haplo.stats
package. This package uses an EM approach that has been adapted to handle estimation of 
haplotype probabilities when the haplotypes are made up of many loci. The haplotype probabilities are estimated for haplotypes
containing all loci. The probability for a given two-locus haplotype is then computed 
by summing up probabilities for the full haplotypes having the given two-locus haplotype (possible 
haplotypes are 00, 01, 10 or 11). These probabilities are computed for all adjacent pairs of loci. 

When using genotype data, it is recommended that a list of likely or known haplotypes and
an initial configuration be provided to sampletrees in order to improve MCMC mixing. These can optionally be initialized using the
output from haplo.em() if HaploList and InitialHaplos are set to TRUE. The list of likely haplotypes will contain those haplotypes that have probability estimated to be above a threshold (set by `tol'). 
The initial haplotype configuration for all individuals is initialized by sampling 
a configuration based on the estimated posterior
probabilities of each haplotype configuration for each individual. 
}
\value{This function writes the estimated haplotype data to the specified files and returns 
\item{args}{An object of class `pars' with the haplotype options set to those 
specified by the call to this function}
}
\references{Burkett KM, McNeney B, Graham J. Sampletrees and Rsampletrees: sampling gene genealogies conditional on SNP genotype data. Bioinformatics. 32:1580-2, 2016}
\author{Kelly Burkett}

\examples{
\dontrun{
datname=paste(path.package("Rsampletrees"),"/extdata/geno_Theta8_Rho8.txt",sep="")
locname=paste(path.package("Rsampletrees"),"/extdata/locations_Theta8_Rho8.txt",sep="")
weightname=paste(path.package("Rsampletrees"),"/extdata/weights-g.txt",sep="")

runpars=newArgs(DataFile=datname, DataType="g", LocationFile=locname, WeightFile=weightname,
	    RunName="Test-g",FocalPoint=10000)

# This will create the three files"EM-hapfreqs", "EM-initial.dat", "EM-known_haplotypes"
runpars=estimateHap(runpars,"EM-hapfreqs",InitialHaploFile="EM-initial.dat",
HaploListFile="EM-known_haplotypes")

}
}
