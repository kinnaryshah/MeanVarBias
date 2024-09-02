library(SpatialExperiment)
library(Rfast2)
library(BiocParallel)
library(here)
library(scran)

# ---------
# load data
# ---------

sim_dir <- "reps_968_50_50per_1000_0.2_to_3_0.5_to_9"
K=1
fn = here("outputs", "simulations", sim_dir, paste0("spe_simulation_",K,".rds"))
spe <- readRDS(file=fn)
spe <- logNormCounts(spe)

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

file = here("outputs","results","simulations",sim_dir,"spe_MoransI.rds")
saveRDS(spe, file = file)

# -----------
# session information
# -----------

sessionInfo()