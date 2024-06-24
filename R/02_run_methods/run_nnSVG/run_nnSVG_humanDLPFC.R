library(SpatialExperiment)
library(nnSVG)
library(scran)
library(spatialLIBD)
library(here)

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_preprocessed.rds")
spe <- readRDS(file=fn)

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

file = here("outputs","results","spe_humanDLPFC_nnSVG.rds"
saveRDS(spe, file = file)


# -----------
# session information
# -----------

sessionInfo()


