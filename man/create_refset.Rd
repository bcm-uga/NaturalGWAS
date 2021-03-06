% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_pheno.R
\name{create_refset}
\alias{create_refset}
\title{Create a set of weakly linked reference SNPs.}
\usage{
create_refset(window = 101, chr.pos)
}
\arguments{
\item{window}{an integer for the spacing of SNPs in the reference set (should be based on LD
information).}

\item{chr.pos}{a p by 2 matrix where p is the number of SNPs (rows) that contains the genetic map for
the genotype matrix. The first column must contain chromosome number (numeric) and the second column
must contain the SNP positions for each chromosome.}
}
\value{
A vector of reference SNP spaced by window size.
}
\description{
This function creates a reference set for causal variants in the simulation based on a LD window.
}
\examples{
library(naturalgwas)

## Load A. thaliana example
data(A.thaliana)
ref.set <- create_refset(A.thaliana$chrpos, window = 101)
}
\seealso{
\code{\link{simu_pheno}}
}
\author{
Olivier François
}
