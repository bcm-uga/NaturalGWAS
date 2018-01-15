#' This function performs a GWAS based on the linear regression method with known
#' confounders. This is an oracle
#' method when the true confounders are provided as arguments.
#' @title Performs a gwas based on linear regression.
#' @author Olivier François
#' @param genotype an n by p matrix of SNP genotypes where n is the number of individuals (rows) and p is the
#' number of SNPs (columns).
#' @param phenotype a vector of n phenotypic traits. If the phenotypes are generated from
#' \code{\link{simu_pheno}} with the same confounder object, then the method is an oracle method.
#' @param confounder an object of class "confounder" created with \code{\link{create_factor}}.
#' @param K number of factors in the adjustment method.
#' @return a vector of p-values and an estimate of genomic inflation for all SNPs in the genotype matrix.
#' @return pvalues  a vector of length p containing all calibrated p-values for phenotype-genotype association tests.
#' @return gif genomic inflation factor.
#' @seealso \code{\link{simu_pheno}} \code{\link{create_factor}}
#' @examples
#' library(naturalgwas)
#'
#' ## Load A. thaliana example
#' data(A.thaliana)
#' env <- get_climate(A.thaliana$coord)
#' confounder <- create_factor(A.thaliana$genotype, K = 10)
#' ref.set <- create_refset(A.thaliana$chrpos, window = 501)
#' sim <- simu_pheno( A.thaliana$genotype, confounder, env, ref.set,
#' ncausal = 8, effect.size = 20, gxe = 1)
#' pv.oracle <- naturalgwas::oracle(sim$phenotype,
#'                                  A.thaliana$genotype,
#'                                  confounder, K = 8)$pv
#'
#' ## compare with 'cate'
#' library(cate)
#' pv.cate = cate::cate( ~ sim.phenotype,
#'                       X.data = data.frame(sim$phenotype),
#'                       Y = A.thaliana$genotype,
#'                       r = 8, calibrate = TRUE)$beta.p.value
#'
#' ## qqplot
#' qqplot(-log10(pv.oracle), -log10(pv.cate) , pch = 19, cex = .4, col = "grey" )
#' abline(0, 1, lwd = 2, col = "orange")
#' par(mfrow = c(2, 1))
#' plot -log10(pv.oracle), cex = .4, col = "grey", main = "Manhattan plot")
#' points( sim$causal, -log10(pv.oracle)[sim$causal], type = "h", lty = 1, col = "blue")
#' plot(-log10(pv.cate), cex = .4, col = "grey", main = "Manhattan plot")
#' points(sim$causal, -log10(pv.cate)[sim$causal], type = "h", lty = 1, col = "red")
#' @export
oracle <- function( phenotype, genotype, confounder, K ){

  if ( K > ncol(confounder$factors) ) {
    stop("K must be lower (or equal) than the dimension of the confounding factors.")
  }

  mod <- lm( as.matrix(genotype) ~ phenotype + confounder$factors[ , 1:K ] )
  z.score <- sapply( summary(mod), function(X) X$coeff[2, 3] )

  gif <- median( z.score^2 ) / qchisq( 0.5, df = 1 )
  p.values <- pchisq( z.score^2 / gif, df = 1, low = FALSE )

  return( list( pvalues = p.values , gif = gif ))
}




