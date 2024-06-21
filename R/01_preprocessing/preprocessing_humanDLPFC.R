###################################
# dataset: Visium human DLPFC
###################################

library(SpatialExperiment)
library(scran)
library(spatialLIBD)
library(nnSVG)
library(here)


# ---------
# load data
# ---------

# load data as a SpatialExperiment from spatialLIBD (already contains logcounts)
spe <- fetch_data(type = "spe")

spe <- spe[, spe$sample_id == "151507"]
spe <- spe[, !is.na(unfactor(colData(spe)$spatialLIBD))]

# -------------
# preprocessing
# -------------

# keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]

dim(spe)


# filter low-expressed 
# using gene filtering function from nnSVG package
spe <- filter_genes(
  spe, 
  filter_genes_ncounts = 2, 
  filter_genes_pcspots = 0.2
)

dim(spe)


# -----------
# save object
# -----------

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_preprocessed.rds")
saveRDS(spe, file = fn)