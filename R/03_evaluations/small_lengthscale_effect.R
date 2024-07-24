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
dat <- c("Ovarian", "DLPFC", "LC", "Breast", "HPC")

# create a table of the output of the function for each dataset
results <- data.frame()
for (d in dat) {
  spe_unweighted <- readRDS(here("outputs", "results", paste0("spe_human", d, "_nnSVG.rds")))
  spe_weighted <- readRDS(here("outputs", "results", paste0("spe_human", d, "_weighted_nnSVG.rds")))
  results <- rbind(results, c(d, small_lengthscale_effect(spe_unweighted, spe_weighted, 0.011, 0.025)))
}

results_old <- results

colnames(results) <- c("Dataset", "Higher Rank", "Total Small Lengthscale", "Prop Unweighted SVGs", "Prop Weighted SVGs", "Total Genes")
results$`Prop Higher Rank` <- results$`Higher Rank` / results$`Total Small Lengthscale`
# move `Prop Higher Rank` to after `Total Small Lengthscale`
results <- results[, c(1, 2, 3, 7, 4, 5, 6)]

# add a column for the proportion of genes that had their weight stabilized
results$`Percent Weight Stabilized` <- c(24.2590865927539, 32.4546892587122, 42.0175677674449, 23.5156050403311, 28.8172894132482)

#   Dataset       Higher Rank     Total Small Lengthscale Prop Higher Rank   Prop Unweighted SVGs Prop Weighted SVGs Total Genes Percent Weight Stabilized
# 1 Ovarian          78                     161             0.4844720           0.3556813          0.3042755       12022                    24.25909
# 2   DLPFC          42                     250             0.1680000           0.2189841          0.2945663        7343                    32.45469
# 3      LC           9                      23             0.3913043           0.5093914          0.5484598        1331                    42.01757
# 4  Breast          77                     106             0.7264151           0.8055353          0.3106891       12321                    23.51561
# 5     HPC          23                      54             0.4259259           0.9280105          0.9820494        5348                    28.81729


#Ovarian 24.2590865927539% of observations had their weight stabilized
#DLPFC 32.4546892587122% of observations had their weight stabilized
#LC 42.0175677674449% of observations had their weight stabilized
#Breast 23.5156050403311% of observations had their weight stabilized
#HPC 28.8172894132482% of observations had their weight stabilized
#Lobular Breast