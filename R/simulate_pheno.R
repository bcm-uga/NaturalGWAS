#' This function downloads the 'wordlclim' bioclimatic database, extract the data for the specified geographic coordinates,
#' and summarizes them with their first principal component.
#' @title Evaluates environment from 'worldclim'
#' @author Olivier François
#' @param coord a matrix of geographic coordinates written as (longitude, latitude).
#' @return the first principal component of 19 bioclimatic variables from 'worldclim'.
#' @seealso \code{\link{simu_pheno}}
#' @examples
#' library(naturalgwas)
#'
#' ## Load example
#' data(A.thaliana)
#' env <- get_climate(A.thaliana$coord)
#' summary(env)
#' @export
get_climate <- function(coord){
  require(raster)
  coord <- as.matrix(coord)
  if (dim(coord)[2] != 2)
    {stop("coord must be a 2 dimensional matrix of longitude and latitude coordinates.")}
  if (dim(coord)[1] == 1)
    {stop("coord must contain more than one point.")}
  climate <-  raster::getData('worldclim', var='bio', res = 2.5)
  bio <-  extract(climate, y = coord)
  pc.bio <- prcomp(bio, scale = T)
  env <- pc.bio$x[,1]
  return(env)
}


#' This function creates confounding factors based on genotypes. Those factors may be included
#' in the simulation of phenotypes.
#' @title Generates confounders.
#' @author Olivier François
#' @param genotype an n by p matrix of genotypes where n is the number of individuals (rows) and p is the
#' number of SNPs (columns).
#' @param K an integer for the number of confounders.
#' @return A \code{list} containing K factors, an estimate of residual error, and an estimate of
#' the order of effects.
#' @return factors a n by K matrix of confounders corresponding to the first K unscaled principal
#' components of the genotype matrix.
#' @return sigma standard deviation for residual errors.
#' @return base order of magnitude for baseline effects.
#' @seealso \code{\link{simu_pheno}}
#' @examples
#' library(naturalgwas)
#'
#' ## Load A. thaliana example
#' data(A.thaliana)
#' confounder <- create_factor(A.thaliana$genotype, K = 20)
#' plot(confounder$factors)
#' cat("Residual error = ", confounder$sigma, "\n")
#' cat("Base effect = ", confounder$base, "\n")
#' @export
create_factor <- function(genotype, K){
  require(LEA)
  LEA::write.geno(R = genotype, output = "genotype.geno")
  pc = LEA::pca("./genotype.geno")
  pc.sdev2 <- pc$sdev^2
  plot(pc.sdev2, lwd = 4, col = "blue", type = "h",
       xlab = "Factors", ylab = "Eigenvalues")

  sigma <- sqrt(sum(pc$sdev^2) - sum(pc$sdev[1:K]^2))
  base.effect <- sqrt(sum(pc$sdev^2))

  result <- list(factors = pc$projections[,1:K],
                 sigma = sigma,
                 base = base.effect)
  remove.pcaProject("genotype.pcaProject")
  file.remove("genotype.geno")
  class(result) <- "confounder"
  return(result)
}


#' This function creates a reference set for causal variants in the simulation based on a LD window.
#' @title Create a set of weakly linked reference SNPs.
#' @author Olivier François
#' @param chr.pos a p by 2 matrix where p is the number of SNPs (rows) that contains the genetic map for
#' the genotype matrix. The first column must contain chromosome number (numeric) and the second column
#' must contain the SNP positions for each chromosome.
#' @param window an integer for the spacing of SNPs in the reference set (should be based on LD
#' information).
#' @return A vector of reference SNP spaced by window size.
#' @seealso \code{\link{simu_pheno}}
#' @examples
#' library(naturalgwas)
#'
#' ## Load A. thaliana example
#' data(A.thaliana)
#' ref.set <- create_refset(A.thaliana$chrpos, window = 101)
#' @export
create_refset <- function(window = 101, chr.pos){
  #definition of a set of reference snps
  lch <- 0
  ref.set <- NULL
  chr <- unique(chr.pos[ , 1])

  for (i in chr){
  chromosome <- which( chr.pos[,1] == i )
  set <-  seq(lch + 1 + (window - 1)/2, lch + length(chromosome), by = window)
  ref.set <- c(ref.set, set)
  lch <- lch + length(chromosome)
  }
  return(ref.set)
}


#' This function samples causal variants from a reference set of variants, and simulates phenotypes based on
#' genotypes, enviromnental variables, and confounding factors.
#' @title Simulate phenotypes based on genotypes.
#' @author Olivier François
#' @param genotype an n by p matrix of SNP genotypes where n is the number of individuals (rows) and p is the
#' number of SNPs (columns).
#' @param ref.set a set of weakly linked reference SNPs created with \code{\link{create_refset}}.
#' @param confounder an object of class "confounder" created with \code{\link{create_factor}}.
#' @param environment an environmental vector created \code{\link{get_climate}}.
#' @param ncausal number of causal variants in the simulation.
#' @param effect.size effect size for causal variants in the simulation.
#' @param gxe intensity of gene by environment interaction in the simulation.
#' @return A subset of 'ncausal' causal variants and a vector of phenotypes.
#' @return causal.set  subset of 'ncausal' causal variants.
#' @return phenotype  vector of phenotypes.
#' @seealso \code{\link{create_refset}} \code{\link{create_factor}} \code{\link{get_climate}}
#' @examples
#' library(naturalgwas)
#'
#' ## Load A. thaliana example
#' data(A.thaliana)
#' env <- get_climate(A.thaliana$coord)
#' confounder <- create_factor(A.thaliana$genotype, K = 20)
#' ref.set <- create_refset(A.thaliana$chrpos, window = 101)
#' sim.ph <- simu_pheno(A.thaliana$genotype,
#'                      confounder,
#'                      env, ref.set,
#'                      ncausal = 8,
#'                      effect.size = 20,
#'                      gxe = 1)
#' @export
simu_pheno <- function(genotype, confounder, environment = NULL, ref.set, ncausal, effect.size, gxe){

  causal.set <- sort( sample(ref.set, ncausal, rep = FALSE) )

  if (!inherits(confounder, "confounder")) stop("The confounder argument is not of class 'confounder'.")

  effect.size <- effect.size*confounder$base
  gxe <- gxe*confounder$base

  if (is.null(environment)){
    if (length(causal.set) == 1){
      pheno <- effect.size*genotype[,causal.set] + rowSums(confounders$factors) + rnorm(nrow(genotype), sd = confounder$sigma)
    }else{
      pheno <- effect.size*rowSums(genotype[ , causal.set]) + rowSums(confounder$factors) + rnorm(nrow(genotype), sd = confounder$sigma)
    }
  } else {
    if (length(causal.set) == 1){
      pheno <- effect.size*genotype[ , causal.set] + gxe*environment*genotype[ , causal.set] + rowSums(confounder$factors) + rnorm(nrow(genotype), sd = confounder$sigma)
    }else{
      pheno <- effect.size*rowSums(genotype[ , causal.set]) + gxe*environment*rowSums(genotype[ , causal.set]) + rowSums(confounder$factors) + rnorm(nrow(genotype), sd = confounder$sigma)
    }
  }

  return(list(phenotype = pheno, causal.set = causal.set))
}
