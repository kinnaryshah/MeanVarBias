library(SpatialExperiment)
library(here)
library(dplyr)

# make a function for this code that takes in as input the spe_unweighted and spe_weighted and low and high lengthscales
# and returns the proportion of genes that have higher ranks after weighting
small_lengthscale_effect <- function(spe_unweighted, spe_weighted, low, high) {
  idx <- which(rowData(spe_unweighted)$phi > low & rowData(spe_unweighted)$phi < high)
  ranks <- data.frame(
    gene = rownames(spe_unweighted)[idx],
    unweighted = rowData(spe_unweighted)$rank[idx],
    weighted = rowData(spe_weighted)$weighted_rank[idx],
    diff = rowData(spe_unweighted)$rank[idx] - rowData(spe_weighted)$weighted_rank[idx]
  )
  print(paste0("prop of unweighted SVGs: ", sum(rowData(spe_unweighted)$padj < 0.05) /dim(spe_unweighted)[1],
               ", prop of weighted SVGs: ", sum(rowData(spe_weighted)$weighted_padj < 0.05, na.rm = T) /dim(spe_weighted)[1],
               ", total genes: ", dim(spe_unweighted)[1]))
  
  return(list(sum(ranks$diff > 0),
              dim(ranks)[1],
              sum(rowData(spe_unweighted)$padj < 0.05) /dim(spe_unweighted)[1],
              sum(rowData(spe_weighted)$weighted_padj < 0.05, na.rm = T) /dim(spe_weighted)[1],
              dim(spe_unweighted)[1]))
}

# list of datasets to use function on
dat <- c("Ovarian", "DLPFC", "LC", "Breast", "HPC", "LobularBreast")

# create a table of the output of the function for each dataset
results <- data.frame()
for (d in dat) {
  spe_unweighted <- readRDS(here("outputs", "results", paste0("spe_human", d, "_nnSVG.rds")))
  spe_weighted <- readRDS(here("outputs", "results", paste0("spe_human", d, "_noise_weighted_nnSVG.rds")))
  results <- rbind(results, c(d, small_lengthscale_effect(spe_unweighted, spe_weighted, 0.011, 0.025)))
}

colnames(results) <- c("Dataset", "Higher Rank", "Total Small Lengthscale", "Prop Unweighted SVGs", "Prop Weighted SVGs", "Total Genes")
results$`Prop Higher Rank` <- results$`Higher Rank` / results$`Total Small Lengthscale`
# move `Prop Higher Rank` to after `Total Small Lengthscale`
results <- results[, c(1, 2, 3, 7, 4, 5, 6)]

# add a column for the proportion of genes that had their weight stabilized
results$`Percent Weight Stabilized` <- c(24.2762214149168, 32.5231071802586, 41.9634057000962, 23.4363301754294, 29.2201329157317, 27.1942350381327)

# Dataset            Higher Rank    Total Small Lengthscale   Prop Higher Rank
# 1       Ovarian          75                     161        0.4658385
# 2         DLPFC          87                     250        0.3480000
# 3            LC           9                      23        0.3913043
# 4        Breast          71                     106        0.6698113
# 5           HPC          21                      54        0.3888889
# 6 LobularBreast          93                     143        0.6503497
# 
# Prop Unweighted SVGs     Prop Weighted SVGs   Total Genes       Percent Weight Stabilized
# 1            0.3556813          0.3785560       12022                  24.27622
# 2            0.2189841          0.2090426        7343                  32.52311
# 3            0.5093914          0.5326822        1331                  41.96341
# 4            0.8055353          0.3808944       12321                  23.43633
# 5            0.9280105          0.9526926        5348                  29.22013
# 6            0.8220057          0.5201996       12624                  27.19424

#Ovarian 24.2762214149168% of observations had their weight stabilized
#DLPFC 32.5231071802586% of observations had their weight stabilized
#LC 41.9634057000962% of observations had their weight stabilized
#Breast 23.4363301754294% of observations had their weight stabilized
#HPC 29.2201329157317% of observations had their weight stabilized
#Lobular Breast 27.1942350381327% of observations had their weight stabilized

