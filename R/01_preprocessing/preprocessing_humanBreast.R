library(scran)
library(SpatialExperiment)
library(here)
# BiocManager::install("DropletUtils") if not installed

(spe <- read10xVisium("outputs/raw/humanBreast/outs",
                      type = "sparse",   # use sparse (not HDF5) format
                      data = "raw",     
                      images = "lowres", # specify which image(s) to include
                      load = TRUE))      # specify whether or not to load image(s)

#row data formatting
rowData(spe)$gene_id <- rownames(spe)
rowData(spe)$gene_name <- rowData(spe)$symbol
rowData(spe)$symbol <- NULL

#remove 94 spots not in tissue
spe <- spe[, colData(spe)$in_tissue == 1]

## remove genes without enough data
spe <- filter_genes(spe, 2, 0.2)

# create logcounts matrix
spe <- logNormCounts(spe)

fn <- here("outputs", "preprocessed", "spe_humanBreast_preprocessed.rds")
saveRDS(spe, file = fn)
