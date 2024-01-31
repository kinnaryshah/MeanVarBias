library(SpatialExperiment)
library(nnSVG)
library(BRISC)
library(BiocParallel)
library(scuttle)
library(spoon)

spe <- spe <- readRDS(file = "spe_simulation.rds")
spe <- spe[, colSums(counts(spe)) > 0]
spe <- logNormCounts(spe)
weights <- generate_weights(input = spe, stabilize = TRUE,
                            BPPARAM = MulticoreParam(workers = 10,
                                                     RNGseed = 4))
save(weights, file = "weights.Rdata")

spe_unweighted <- nnSVG(spe, assay_name = "logcounts", n_threads = 10)
spe_weighted <- weighted_nnSVG(input = spe, w = weights,
                            BPPARAM = MulticoreParam(workers = 10,
                                                     RNGseed = 6))

save(spe_unweighted, spe_weighted, file = "spe_simulation_weighted_nnSVG.Rdata")
