## ------------------------------------------------------------------------
#install.packages("cate")
#devtools::install_github("bcm-uga/lfmm")
library(naturalgwas)

## ----data----------------------------------------------------------------
data(A.thaliana)
genotype <- A.thaliana$genotype
coordinates <- A.thaliana$coord
chrpos <- A.thaliana$chrpos 

## ------------------------------------------------------------------------
dim(genotype)
genotype[1:3,1:10]

## ------------------------------------------------------------------------
library(maps)
plot(coordinates, pch = 19, cex = .5, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
map(add = T, interior = F)

## ----init, dependson=c("data"), results="hide"---------------------------
  env <- get_climate( coordinates )
  ref.set <- create_refset( chrpos, window = 101 )
  confounder <- create_factor( genotype, K = 20 )

## ---- dependson=c("init")------------------------------------------------
  confounder$sigma 
  confounder$base 

## ----sim, dependson=c("init"), results="hide"----------------------------
  sim <- simu_pheno( A.thaliana$genotype, confounder, env, ref.set,
                     ncausal = 20, effect.size = 500, gxe = .3 )

## ---- dependson=c("sim")-------------------------------------------------
  hist(sim$phenotype, main = "Trait value")

## ---- dependson=c("sim")-------------------------------------------------
  library(fields)
  fit = fields::Krig(coordinates, sim$phenotype, theta = 10, m = 2)
  surface(fit, extrap = TRUE, xlab = "Longitude", ylab = "Latitude", levels = c(0))
  map(add = TRUE, interior = F)

## ----oracle, dependson=c("sim")------------------------------------------
pv.oracle <- oracle( sim$phenotype, genotype, confounder, K = 6 )$pv

## ----cate, dependson=c("oracle")-----------------------------------------
  library(cate)
  pv.cate <- cate( ~ sim.phenotype, X.data = data.frame(sim$phenotype), Y = genotype, r = 6, calibrate = TRUE)$beta.p.value

## ---- dependson=c("cate")------------------------------------------------
   plot( -log10(pv.oracle), cex = .4, col = "grey", main = "Manhattan plot (Oracle)" )
   points( sim$causal, -log10(pv.oracle)[sim$causal], type = "h", lty = 1, col = "blue" )
   plot( -log10(pv.cate), cex = .4, col = "grey", main = "Manhattan plot (Latent factor model)" )
   points( sim$causal, -log10(pv.cate)[sim$causal], type = "h", lty = 1, col = "orange" )

## ---- dependson=c("cate")------------------------------------------------
   ## qqplot
   qqplot( -log10(pv.oracle), -log10(pv.cate) , pch = 19, cex = .4, col = "grey" )
   abline( 0, 1, lwd = 2, col = "orange" )

## ----lfmm----------------------------------------------------------------
lfmm.res <- lfmm::lfmm_ridge(Y = scale(genotype, scale = FALSE), X = scale(sim$phenotype), K = 6)

## ------------------------------------------------------------------------
  p <- lfmm::lfmm_test(Y = scale(genotype, scale = FALSE), X = scale(sim$phenotype), lfmm = lfmm.res, calibrate = "gif")
  pv.lfmm <- p$calibrated.pvalue

## ---- dependson=c("cate")------------------------------------------------
   plot( -log10(pv.lfmm), cex = .4, col = "grey", main = "Manhattan plot (lfmm)" )
   points( sim$causal, -log10(pv.lfmm)[sim$causal], type = "h", lty = 1, col = "blue" )
   plot( -log10(pv.cate), cex = .4, col = "grey", main = "Manhattan plot (cate)" )
   points( sim$causal, -log10(pv.cate)[sim$causal], type = "h", lty = 1, col = "orange" )

## ---- dependson=c("cate")------------------------------------------------
   ## plot
   plot(-log10(pv.lfmm), -log10(pv.cate) , pch = 19, cex = .4, col = "grey" )
   abline( 0, 1, lwd = 2, col = "blue" )

