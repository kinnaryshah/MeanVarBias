library(here)
library(SpatialExperiment)
library(nnSVG)
library(scran)

# downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211956
# created similar SpaceRanger-like directory structure

spe <- read10xVisium("~/ovarian_cancer_forrest/outs")

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
