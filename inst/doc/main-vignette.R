## ------------------------------------------------------------------------
#install.packages("cate")
#devtools::install_github("bcm-uga/lfmm")
library(naturalgwas)

## ----data----------------------------------------------------------------
data(A.thaliana)
genotype <- A.thaliana$genotype
chrpos <- A.thaliana$chrpos 

## ------------------------------------------------------------------------
dim(genotype)
genotype[sample(162, 3),1:10]

## ------------------------------------------------------------------------
library(maps)
coordinates <- A.thaliana$coord
plot(coordinates, 
     pch = 19, cex = .7, col = "green4",
     xlab = "Longitude (°E)", 
     ylab = "Latitude (°N)")
map(add = T, interior = FALSE)

## ----init, results="hide"------------------------------------------------
  env <- get_climate(coordinates)
  ref.set <- create_refset(chrpos, window = 501)
  confounder <- create_factor(genotype, K = 10)

## ---- dependson=c("init")------------------------------------------------
  confounder$sigma 
  confounder$base 

## ----sim, dependson=c("init"), results="hide"----------------------------
  n.causal <- 6
  eff.size <- 5.0

  sim <- simu_pheno(A.thaliana$genotype, 
                    confounder, 
                    env, 
                    ref.set, ncausal = n.causal, 
                    effect.size = eff.size, gxe = .5)
  
  ss <- sum(apply(genotype[,sim$causal.set], 2, var))
  ss*(eff.size)^2*confounder$base^2/(ss*(eff.size)^2*confounder$base^2 + confounder$sigma^2 + confounder$base^2)
  (ss*(eff.size)^2*confounder$base^2 + confounder$sigma^2)/(ss*(eff.size)^2*confounder$base^2 + confounder$sigma^2 + confounder$base^2)

## ---- dependson=c("sim")-------------------------------------------------
  hist(scale(sim$phenotype), main = "Trait value")

## ---- dependson=c("sim")-------------------------------------------------
  library(fields)
  fit = fields::Krig(coordinates, scale(sim$phenotype), theta = 10, m = 2)
  surface(fit, extrap = TRUE, xlab = "Longitude", ylab = "Latitude", levels = c(0))
  map(add = TRUE, interior = F)

## ----oracle, dependson=c("sim")------------------------------------------
pv.oracle <- oracle(scale(sim$phenotype), 
                    genotype, 
                    confounder, K = 10)$pv

## ---- dependson=c("cate")------------------------------------------------
   plot(-log10(pv.oracle), cex = .4, col = "grey", main = "Manhattan plot (Oracle)")
   points(sim$causal, -log10(pv.oracle)[sim$causal], type = "h", lty = 1, col = "blue")
   abline(-log10(0.05/length(pv.oracle)), 0, lty = 2)

## ------------------------------------------------------------------------
library(rrBLUP)
G <- A.thaliana$genotype
rownames(G) <-  A.thaliana$ecotype.id[,1]
data <- data.frame(y = scale(sim$phenotype), 
                   gid = A.thaliana$ecotype.id[,1])

#predict breeding values and compute h2
kb <- kin.blup(data = data, 
                geno="gid", 
                pheno = "y", 
                K = A.mat(G))

cat("Heritability = ", kb$Vg/(kb$Vg + kb$Ve), "\n")

## ------------------------------------------------------------------------
## fit latent factors using an LFMM
mod.lfmm <- lfmm::lfmm_ridge(Y = A.thaliana$genotype, 
                             X = scale(sim$phenotype), 
                             K = 6)

Kinship <- mod.lfmm$U %*% t(mod.lfmm$U)/nrow(mod.lfmm$U)
rownames(Kinship) <- A.thaliana$ecotype.id[,1]
## predict breeding values and compute h2
kb <- kin.blup(data = data, 
                geno="gid", 
                pheno = "y", 
                K = Kinship)

cat("Heritability = ", kb$Vg/(kb$Vg + kb$Ve), "\n")

## ------------------------------------------------------------------------
geno <- t(A.thaliana$genotype)
colnames(geno) <- A.thaliana$ecotype.id[,1]
G <- data.frame(marker = 1:nrow(A.thaliana$chrpos),
                A.thaliana$chrpos,
                geno,
                check.names = FALSE)
pheno <- data.frame(line = A.thaliana$ecotype.id[,1], 
                    y = scale(sim$phenotype))

mlm <- GWAS(pheno = pheno, 
            geno = G, 
            K = Kinship,
            n.PC = 5,    
            plot = FALSE)

## ------------------------------------------------------------------------
## Manhattan plot
plot(mlm$y, cex = .4, col = "grey",
     ylab = "-log10(pvalues)",
     main = "Manhattan plot (mlm)" )
points(sim$causal, mlm$y[sim$causal], 
       type = "h", lty = 1, col = "green4")
abline(-log10(0.05/length(pv.oracle)), 0, lty = 2)

## ------------------------------------------------------------------------
qqplot(-log10(pv.oracle), mlm$y, col = "grey", cex = .4, ylab ="mlm scores")
abline(0,1, col = 'green4')

## ----cate, dependson=c("oracle")-----------------------------------------
  pheno <- scale(sim$phenotype)
  pv.cate <- cate::cate( ~ pheno, 
                   X.data = data.frame(pheno), 
                   Y = genotype, 
                   r = 6, 
                   calibrate = TRUE)$beta.p.value

## ---- dependson=c("cate")------------------------------------------------
   plot(-log10(pv.cate), cex = .4, col = "grey", main = "Manhattan plot (Latent factor model)" )
   points(sim$causal, -log10(pv.cate)[sim$causal], type = "h", lty = 1, col = "orange" )
    abline(-log10(0.05/length(pv.oracle)), 0, lty = 2)

## ---- dependson=c("cate")------------------------------------------------
   ## qqplot
   qqplot(-log10(pv.oracle), -log10(pv.cate) , pch = 19, cex = .4, col = "grey" )
   abline( 0, 1, lwd = 2, col = "orange" )

## ---- dependson=c("cate")------------------------------------------------
   ## qqplot
   qqplot(mlm$y, -log10(pv.cate) , pch = 19, cex = .4, col = "grey" )
   abline( 0, 1, lwd = 2, col = "orange" )

## ----lfmm----------------------------------------------------------------
lfmm.mod <- lfmm::lfmm_ridge(Y = genotype, 
                             X = scale(sim$phenotype), 
                             K = 6)

## ------------------------------------------------------------------------
  p <- lfmm::lfmm_test(Y = genotype, 
                       X = scale(sim$phenotype), 
                       lfmm = lfmm.mod, 
                       calibrate = "gif")
  pv.lfmm <- p$calibrated.pvalue

## ---- dependson=c("cate")------------------------------------------------
   plot( -log10(pv.lfmm), cex = .4, col = "grey", main = "Manhattan plot (lfmm)" )
   points( sim$causal, -log10(pv.lfmm)[sim$causal], type = "h", lty = 1, col = "blue" )
   abline(-log10(0.05/length(pv.oracle)), 0, lty = 2)

## ---- dependson=c("cate")------------------------------------------------
   ## plot
   qqplot(-log10(pv.oracle), -log10(pv.lfmm) , pch = 19, cex = .4, col = "grey" )
   abline( 0, 1, lwd = 2, col = "blue" )

## ------------------------------------------------------------------------
cond.test <- lfmm::forward_test(Y = A.thaliana$genotype,
                                X = scale(sim$phenotype),
                                K = 6,
                                niter = 12)

## ------------------------------------------------------------------------
# proposed candidates
cat("Candidate loci:")
sort(cond.test$candidates)
# truth
cat("Causal SNPs:")
sim$causal.set

