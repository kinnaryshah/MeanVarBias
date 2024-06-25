library(SpatialExperiment)
library(nnSVG)
library(scran)
library(spatialLIBD)
library(here)

# ---------
# load data
# ---------

fn <- here("outputs", "preprocessed", "spe_humanBreast_preprocessed.rds")
spe <- readRDS(file=fn)

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

file = here("outputs","results","spe_humanBreast_nnSVG.rds")
saveRDS(spe, file = file)


# -----------
# session information
# -----------

sessionInfo()


