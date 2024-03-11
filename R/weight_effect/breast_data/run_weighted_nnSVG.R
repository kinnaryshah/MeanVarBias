library(nnSVG)
library(SpatialExperiment)
library(BRISC)
library(BiocParallel)
library(scuttle)
library(Matrix)
library(spoon)


(spe <- read10xVisium("/users/kshah/mean_var_project/breast_cancer_data/outs",
                      type = "sparse",   # use sparse (not HDF5) format
                      data = "raw",     
                      images = "lowres", # specify which image(s) to include
                      load = TRUE))      # specify whether or not to load image(s)

#row data did not load in the proper format
rowData(spe)$gene_id <- rownames(spe)
rowData(spe)$gene_name <- rowData(spe)$symbol
rowData(spe)$symbol <- NULL

dim(spe)

#removes 94 spots not in tissue
spe <- spe[, colData(spe)$in_tissue == 1]

dim(spe)

#convert counts matrix from sparse matrix to matrix
#weighted_nnSVG requires a matrix for calculating weighted_logcounts matrix
counts(spe) <- as.matrix(counts(spe))

## Remove genes without enough data
spe <- filter_genes(spe, 2, 0.2)

dim(spe)

spe <- logNormCounts(spe)

weights <- generate_weights(input = spe,
                            stabilize = TRUE,
                            BPPARAM = MulticoreParam(workers = 10,
                                                     RNGseed = 4))
print(dim(weights))
saveRDS(weights, file = "weights.rds")

spe <- weighted_nnSVG(input = spe,
                      w = weights,
                      BPPARAM = MulticoreParam(workers = 10, RNGseed = 5))

file = "spe_weighted_nnSVG.rds"
saveRDS(spe, file = file)
