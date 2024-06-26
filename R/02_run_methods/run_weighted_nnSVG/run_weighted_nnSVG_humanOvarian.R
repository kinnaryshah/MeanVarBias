library(SpatialExperiment)
library(nnSVG)
library(scran)
library(here)
library(spoon)
library(Matrix)
library(scuttle)
library(BRISC)
library(BiocParallel)

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
# run weighted nnSVG
# ---------

weights <- generate_weights(
  input = spe,
  stabilize = TRUE,
  BPPARAM = MulticoreParam(workers = 10, RNGseed = 4))

spe <- weighted_nnSVG(
  input = spe,
  w = weights,
  BPPARAM = MulticoreParam(workers = 10, RNGseed = 5))

# -----------
# save objects
# -----------

file = here("outputs","results","spe_humanOvarian_weighted_nnSVG.RData")
save(spe, weights, file = file)

# -----------
# session information
# -----------

sessionInfo()


