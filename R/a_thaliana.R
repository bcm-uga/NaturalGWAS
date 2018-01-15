#' Genetic and geographic data for Arabidopsis thaliana
#'
#' A dataset containing SNP frequency and geographic data for 162 European
#' plant accessions from Atwell et al. Science 2010.
#' The variables are as follows:
#'
#' \itemize{
#'   \item genotype: binary (0 or 1) SNP frequency for 162 European accessions (53859 SNPs).
#'   \item coordinates: longitude and latitude of individuals.
#'   \item chrpos: genetic map, chromosome (Chr 5) and position for each SNP.
#'   \item f16: days to flowering trait (FT16) from Atwell et al. 2010. Plants
#'   were checked bi-weekly for presence of first buds, and the average flowering time of
#'   4 plants of the same accession were collected.
#'  \item ecotype.id: ecotype.id, native name, site and region for 162 European accessions.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name A.thaliana
#' @usage data(A.thaliana)
#' @format A list with 4 arguments: genotype, coord, chrpos, f16
NULL
