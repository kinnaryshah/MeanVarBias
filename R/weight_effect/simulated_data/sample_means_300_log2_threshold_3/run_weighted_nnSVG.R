library(nnSVG)
library(SpatialExperiment)
library(BRISC)
library(BiocParallel)
library(scuttle)
library(Matrix)
library(spoon)

fn <- ("../../../simulations/sample_means_300_log2_threshold_3/spe_simulation.rds") 

spe <- readRDS(fn)

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
