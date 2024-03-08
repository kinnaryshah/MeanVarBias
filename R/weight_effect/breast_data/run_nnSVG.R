library(nnSVG)
library(spatialLIBD)
library(scran)

(spe <- read10xVisium("/users/kshah/mean_var_project/breast_cancer_data/outs",
                      type = "sparse",   # use sparse (not HDF5) format
                      data = "raw",     
                      images = "lowres", # specify which image(s) to include
                      load = TRUE))      # specify whether or not to load image(s)

#row data did not load in the proper format
rowData(spe)$gene_id <- rownames(spe)
rowData(spe)$gene_name <- rowData(spe)$symbol
rowData(spe)$symbol <- NULL

#removes 94 spots not in tissue
spe <- spe[, colData(spe)$in_tissue == 1]

dim(spe)

## Remove genes without enough data
spe <- filter_genes(spe, 2, 0.2)

dim(spe)

spe <- logNormCounts(spe)

runtime <- system.time({
  tryCatch({
    spe <- nnSVG(
      spe,
      X = NULL,
      assay_name = "logcounts",
      n_neighbors = 10,
      order = "AMMD",
      n_threads = 10,
      verbose = FALSE
    )
  }, error = function(e) {
    print(paste("Error:", e$message))
  })
})

file = "spe_nnSVG.rds"
saveRDS(spe, file = file)
