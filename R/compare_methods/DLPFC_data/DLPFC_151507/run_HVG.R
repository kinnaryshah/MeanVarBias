library(SpatialExperiment)
library(scran)
library(spatialLIBD)
library(nnSVG)

#already contains logcounts
spe <- fetch_data(type = "spe")

spe <- spe[, spe$sample_id == "151507"]
spe <- spe[, !is.na(unfactor(colData(spe)$spatialLIBD))]

## Remove genes without enough data
spe <- filter_genes(spe, 2, 0.2)

dim(spe)

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
