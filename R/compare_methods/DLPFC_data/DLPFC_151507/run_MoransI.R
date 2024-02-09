library(SpatialExperiment)
library(Rfast2)
library(BiocParallel)
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

# note that this needs logcounts
spe <- logNormCounts(spe)

d <- dist(spatialCoords(spe))
d_mat <- as.matrix(d)
w <- 1 / (d_mat ^ 2)
diag(w) <- 0

# run Moran's I
set.seed(123)
runtime <- system.time({
  res <- bplapply(seq_len(nrow(spe)), function(i) {
    y_i <- logcounts(spe)[i, ]
    out_i <- moranI(
      x = y_i, 
      w = w, 
      scaled = FALSE, 
      R = 0
    )
    unname(out_i[1])
  }, BPPARAM = MulticoreParam(workers = 10))
})

res <- unlist(res)

# store in object
stopifnot(length(res) == nrow(spe))
rowData(spe)$MoransI <- res

# calculate ranks
rowData(spe)$MoransI_rank <- rank(-1 * rowData(spe)$MoransI, ties.method = "first")

file <- "spe_MoransI.rds"
saveRDS(spe, file=file)
