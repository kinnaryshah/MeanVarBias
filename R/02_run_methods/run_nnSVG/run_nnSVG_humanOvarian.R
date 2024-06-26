library(SpatialExperiment)
library(nnSVG)
library(scran)
library(spatialLIBD)
library(here)

# ---------
# load data
# ---------

fn <- here("outputs", "preprocessed", "spe_humanOvarian_preprocessed.rds")
spe <- readRDS(file=fn)

# ---------
# convert to dgCMatrix
# ---------

regular_matrix_counts <- as.matrix(assays(spe)[["counts"]])
sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
assays(spe)$counts <- sparse_matrix_counts

regular_matrix_logcounts <- as.matrix(assays(spe)[["logcounts"]])
sparse_matrix_logcounts <- as(regular_matrix_logcounts, "dgCMatrix")
assays(spe)$logcounts <- sparse_matrix_logcounts

# ---------
# run nnSVG
# ---------

set.seed(123)
spe <- nnSVG(
  spe,
  X = NULL,
  assay_name = "logcounts",
  n_neighbors = 10,
  order = "AMMD",
  n_threads = 10,
  verbose = FALSE
)

# -----------
# save object
# -----------

file = here("outputs","results","spe_humanOvarian_nnSVG.rds")
saveRDS(spe, file = file)


# -----------
# session information
# -----------

sessionInfo()


