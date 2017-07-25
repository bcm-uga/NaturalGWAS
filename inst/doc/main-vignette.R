## ------------------------------------------------------------------------
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
                     ncausal = 20, effect.size = 5000, gxe = .3 )

## ---- dependson=c("sim")-------------------------------------------------
  hist(sim$phenotype, main = "Trait value")

## ---- dependson=c("sim")-------------------------------------------------
  library(fields)
  fit = Krig(coordinates, sim$phenotype, theta = 10, m = 2)
  surface(fit, extrap = TRUE, xlab = "Longitude", ylab = "Latitude")
  map(add = TRUE, interior = F)

## ----oracle, dependson=c("sim")------------------------------------------
pv.oracle <- oracle( sim$phenotype, genotype, confounder, K = 6 )$pv

## ----cate, dependson=c("oracle")-----------------------------------------
  library(cate)
  pv.cate <- cate( ~ sim.phenotype, X.data = data.frame(sim$phenotype), Y = genotype, r = 6, calibrate = TRUE)$beta.p.value

## ---- dependson=c("cate")------------------------------------------------
   plot( -log10(pv.oracle), cex = .4, col = "grey", main = "Manhattan plot" )
   points( sim$causal, -log10(pv.oracle)[sim$causal], type = "h", lty = 1, col = "blue" )
   plot( -log10(pv.cate), cex = .4, col = "grey", main = "Manhattan plot" )
   points( sim$causal, -log10(pv.cate)[sim$causal], type = "h", lty = 1, col = "red" )

## ---- dependson=c("cate")------------------------------------------------
   ## qqplot
   qqplot( -log10(pv.oracle), -log10(pv.cate) , pch = 19, cex = .4, col = "grey" )
   abline( 0, 1, lwd = 2, col = "orange" )

