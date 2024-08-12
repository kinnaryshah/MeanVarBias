###################################
# dataset: Visium human hippocampus
###################################
library(here)
library(SpatialExperiment)
library(nnSVG)
library(scran)
# BiocManager::install("DropletUtils") if not installed

# ---------
# load data
# ---------

spe <- read10xVisium(here("outputs/raw/humanHPC/outs/"),
                     type = "sparse",   # use sparse (not HDF5) format
                     data = "raw",     
                     images = "lowres", # specify which image(s) to include
                     load = TRUE)      # specify whether or not to load image(s)

# -------------
# preprocessing
# -------------

#row data formatting
rowData(spe)$gene_id <- rownames(spe)
rowData(spe)$gene_name <- rowData(spe)$symbol
rowData(spe)$symbol <- NULL

spe <- filter_genes(spe)

spe <- logNormCounts(spe)


# -----------
# save object
# -----------

fn <- here("outputs", "preprocessed", "spe_humanHPC_preprocessed.rds")
saveRDS(spe, file = fn)

# -----------
# session information
# -----------

sessionInfo()

