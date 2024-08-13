library(SpatialExperiment)
library(dplyr)
library(here)
library(ggplot2)
library(escheR)
library(PRECAST)
library(Seurat)
library(tidyverse)

# Lobular Breast cancer related genes
spe_unweighted <- readRDS(here("outputs", "results", "spe_humanLobularBreast_nnSVG.rds"))
genes <- c("PDK1", "PPAT", "ATP6AP1L", "LYRM7", "TRIM35", "TRIM32", "FRAT1", "TIGAR", "ZNF439", "MOSPD2", "BRWD3")
indices <- which(rowData(spe_unweighted)$gene_name %in% genes)

# in colData(spe_unweighted), rename array_row and array_col to row and col
colData(spe_unweighted)$row <- colData(spe_unweighted)$array_row
colData(spe_unweighted)$col <- colData(spe_unweighted)$array_col

# use clustering with k=6 to cluster the cells
seu <- CreateSeuratObject(counts = as.matrix(counts(spe_unweighted)), 
                          meta.data=data.frame(colData(spe_unweighted)),
                          project="LobularBreast")

seuList <- list("LobularBreast" = seu)

set.seed(1)
preobj <- CreatePRECASTObject(seuList = seuList, gene.number=2000, selectGenesMethod='HVGs',
                              premin.spots = 1, premin.features=1, postmin.spots=1, postmin.features=1)
preobj@seulist

PRECASTObj <- AddAdjList(preobj, platform = "Visium")

PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE,  maxIter = 30, verbose = TRUE)

K <- 6
PRECASTObj <- PRECAST(PRECASTObj, K = K)

PRECASTObj <- SelectModel(PRECASTObj)
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

# Merge with spe object
cluster_df <- seuInt@meta.data |>
  mutate(cluster = factor(cluster)) |>
  rename_with(~ paste0("PRECAST_", .x)) |>
  rownames_to_column(var = "key")

colData(spe_unweighted)$key <- rownames(colData(spe_unweighted))

col_data_df <- colData(spe_unweighted) |>
  data.frame() |>
  left_join(cluster_df, by="key")

rownames(col_data_df) <- colnames(spe_unweighted)
colData(spe_unweighted)$PRECAST_cluster <- col_data_df$PRECAST_cluster

genes <- c("PDK1", "PPAT", "ATP6AP1L", "LYRM7", "TRIM35", "TRIM32", "FRAT1", "TIGAR", "ZNF439", "MOSPD2", "BRWD3")

pdf(here("plots", "transcription_factors_LobularBreast_escheR_6.pdf"))
for (gene in genes) {
  spe_unweighted$counts_gene <- counts(spe_unweighted)[which(rowData(spe_unweighted)$gene_name==gene),]
  p <- make_escheR(spe_unweighted)  |> 
    add_fill(var = "counts_gene") |> 
    add_ground(var = "PRECAST_cluster", stroke = 0.5) + 
    scale_fill_gradient(low = "white", high = "black")
  print(p)
}

dev.off()
