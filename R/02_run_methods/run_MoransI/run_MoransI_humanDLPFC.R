library(SpatialExperiment)
library(Rfast2)
library(BiocParallel)
library(here)

# ---------
# load data
# ---------

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_preprocessed.rds")
spe <- readRDS(file=fn)

# ---------
# run MoransI
# ---------

d <- dist(spatialCoords(spe))
d_mat <- as.matrix(d)
w <- 1 / (d_mat ^ 2)
diag(w) <- 0

set.seed(123)
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


res <- unlist(res)

# store results in spe object
stopifnot(length(res) == nrow(spe))
rowData(spe)$MoransI <- res

# calculate ranks
rowData(spe)$MoransI_rank <- rank(-1 * rowData(spe)$MoransI, ties.method = "first")

# -----------
# save object
# -----------

file = here("outputs","results","spe_humanDLPFC_MoransI.rds")
saveRDS(spe, file = file)

# -----------
# session information
# -----------

sessionInfo()