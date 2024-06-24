library(here)
library(SpatialExperiment)
library(nnSVG)
library(scran)
# BiocManager::install("DropletUtils") if not installed

(spe <- read10xVisium("outputs/raw/humanOvarian/outs",
                     type = "sparse",   # use sparse (not HDF5) format
                     data = "filtered",     
                     images = "lowres", # specify which image(s) to include
                     load = TRUE))      # specify whether or not to load image(s)


#row data formatting
rowData(spe)$gene_id <- rownames(spe)
rowData(spe)$gene_name <- rowData(spe)$symbol
rowData(spe)$symbol <- NULL

## Remove genes without enough data
spe <- filter_genes(spe, 2, 0.2)

dim(spe)

spe <- logNormCounts(spe)

fn <- here("outputs", "preprocessed", "spe_humanOvarian_preprocessed.rds")
saveRDS(spe, file = fn)
