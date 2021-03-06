% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwas_methods.R
\name{oracle}
\alias{oracle}
\title{Performs a gwas based on linear regression.}
\usage{
oracle(phenotype, genotype, confounder, K)
}
\arguments{
\item{phenotype}{a vector of n phenotypic traits. If the phenotypes are generated from
\code{\link{simu_pheno}} with the same confounder object, then the method is an oracle method.}

\item{genotype}{an n by p matrix of SNP genotypes where n is the number of individuals (rows) and p is the
number of SNPs (columns).}

\item{confounder}{an object of class "confounder" created with \code{\link{create_factor}}.}

\item{K}{number of factors in the adjustment method.}
}
\value{
a vector of p-values and an estimate of genomic inflation for all SNPs in the genotype matrix.

pvalues  a vector of length p containing all calibrated p-values for phenotype-genotype association tests.

gif genomic inflation factor.
}
\description{
This function performs a GWAS based on the linear regression method with known
confounders. This is an oracle
method when the true confounders are provided as arguments.
}
\examples{
library(naturalgwas)

## Load A. thaliana example
data(A.thaliana)
env <- get_climate(A.thaliana$coord)
confounder <- create_factor(A.thaliana$genotype, K = 10)
ref.set <- create_refset(A.thaliana$chrpos, window = 501)
sim <- simu_pheno( A.thaliana$genotype, confounder, env, ref.set,
ncausal = 8, effect.size = 20, gxe = 1)
pv.oracle <- naturalgwas::oracle(sim$phenotype,
                                 A.thaliana$genotype,
                                 confounder, K = 8)$pv

## compare with 'cate'
library(cate)
pv.cate = cate::cate( ~ sim.phenotype,
                      X.data = data.frame(sim$phenotype),
                      Y = A.thaliana$genotype,
                      r = 8, calibrate = TRUE)$beta.p.value

## qqplot
qqplot(-log10(pv.oracle), -log10(pv.cate) , pch = 19, cex = .4, col = "grey" )
abline(0, 1, lwd = 2, col = "orange")
par(mfrow = c(2, 1))
plot -log10(pv.oracle), cex = .4, col = "grey", main = "Manhattan plot")
points( sim$causal, -log10(pv.oracle)[sim$causal], type = "h", lty = 1, col = "blue")
plot(-log10(pv.cate), cex = .4, col = "grey", main = "Manhattan plot")
points(sim$causal, -log10(pv.cate)[sim$causal], type = "h", lty = 1, col = "red")
}
\seealso{
\code{\link{simu_pheno}} \code{\link{create_factor}}
}
\author{
Olivier François
}
