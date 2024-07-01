library(spatialLIBD)
library(nnSVG)
library(here)

# ---------
# load data
# ---------

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_preprocessed.rds")
spe <- readRDS(file=fn)

# ---------
# run nnSVG for each layer
# ---------

layer_list <- c("L1", "L2", "L3", "L4", "L5", "L6", "WM")

for (layer in layer_list) {
  spe_layer <- spe[, colData(spe)$spatialLIBD == layer]
  
  # filter out genes with low expression
  spe_layer <- filter_genes(spe_layer, 2, 0.2)

  # set seed for reproducibility
  set.seed(123)
  
  #run nnSVG
  spe_layer <- nnSVG(
    spe_layer,
    X = NULL,
    assay_name = "logcounts",
    n_neighbors = 10,
    order = "AMMD",
    n_threads = 10,
    verbose = FALSE
  )
  
  # save results
  file = here("outputs","results",paste0("spe_humanDLPFC_layer_", layer, "_nnSVG.rds"))
  saveRDS(spe_layer, file = file)
}
