library(SpatialExperiment)
library(scran)

fn <- ("../../../simulations/sample_means_300/spe_simulation.rds") 

spe <- readRDS(fn)

dim(spe)

# note that this needs logcounts
spe <- logNormCounts(spe)

# run HVGs
set.seed(123)
runtime <- system.time({
  dec <- modelGeneVar(spe)
})

# store in object
stopifnot(all(rownames(dec) == rowData(spe)$gene_id))
rowData(spe) <- cbind(rowData(spe), dec)

# calculate ranks
rowData(spe)$HVG_rank <- rank(-1 * rowData(spe)$bio, ties.method = "first")

file = "spe_HVG.rds"
saveRDS(spe, file=file)
