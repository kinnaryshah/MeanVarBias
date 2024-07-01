
# devtools::install_github('xzhoulab/SPARK')
library(SpatialExperiment)
library(SPARK)
library(here)

# ---------
# load data
# ---------

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_preprocessed.rds")
spe <- readRDS(file=fn)

# ---------
# run SPARKX
# ---------

set.seed(123)
sparkx_out <- sparkx(
  count_in = counts(spe),
  locus_in = spatialCoords(spe),
  X_in = NULL,
  numCores = 1,
  option = "mixture",
  verbose = TRUE
)


# store results in spe object
stopifnot(all(rownames(sparkx_out$res_mtest) == rowData(spe)$gene_id))

rowData(spe) <- cbind(rowData(spe), sparkx_out$res_mtest)

# calculate ranks
rowData(spe)$SPARKX_rank <- rank(rowData(spe)$combinedPval, ties.method = "first")

# -----------
# save object
# -----------

file = here("outputs","results","spe_humanDLPFC_SPARKX.rds")
saveRDS(spe, file = file)

# -----------
# session information
# -----------

sessionInfo()