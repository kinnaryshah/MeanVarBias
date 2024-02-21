library(nnSVG)
library(SpatialExperiment)
library(BRISC)
library(BiocParallel)
library(scuttle)
library(Matrix)
library(spoon)

fn <- (file = "mean_var_project/spe_V10B01_086_B1_HPC.rds")

spe <- readRDS(fn)

dim(spe)

## Remove genes without enough data
spe <- filter_genes(spe, 2, 0.2)

dim(spe)

spe <- logNormCounts(spe)

weights <- generate_weights(input = spe,
                            stabilize = TRUE,
                            BPPARAM = MulticoreParam(workers = 10,
                                                     RNGseed = 4))

spe <- weighted_nnSVG(input = spe,
                      w = weights,
                      BPPARAM = MulticoreParam(workers = 10, RNGseed = 5))

file = "spe_weighted_nnSVG.rds"
saveRDS(spe, file = file)
