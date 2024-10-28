library(SpatialExperiment)
library(here)
library(dplyr)
library(ggplot2)

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
dat <- c("Ovarian", "DLPFC", "LC", "Breast", "HPC_V12D07-335_D1", "LobularBreast", "SubtypeBreast")

# create a table of the output of the function for each dataset
results <- data.frame()
for (d in dat) {
  spe_unweighted <- readRDS(here("outputs", "results", paste0("spe_human", d, "_nnSVG.rds")))
  spe_weighted <- readRDS(here("outputs", "results", paste0("spe_human", d, "_weighted_nnSVG.rds")))
  results <- rbind(results, c(d, small_lengthscale_effect(spe_unweighted, spe_weighted, 0.011, 0.025)))
}

colnames(results) <- c("Dataset", "Higher Rank", "Total Small Lengthscale", "Prop Unweighted SVGs", "Prop Weighted SVGs", "Total Genes")
results$`Prop Higher Rank` <- results$`Higher Rank` / results$`Total Small Lengthscale`

# create a plot with x-axis as prop unweighted SVGs and y-axis as prop higher rank
# rename HPC_V12D07-335_D1 to HPC
results$Dataset[results$Dataset == "HPC_V12D07-335_D1"] <- "HPC"
# rename Breast to Ductal Breast
results$Dataset[results$Dataset == "Breast"] <- "Ductal Breast"
# rename SubtypeBreast to Subtype Breast
results$Dataset[results$Dataset == "SubtypeBreast"] <- "Subtype Breast"
# rename LobularBreast to Lobular Breast
results$Dataset[results$Dataset == "LobularBreast"] <- "Lobular Breast"

fig1 <- ggplot(results, aes(x = `Prop Unweighted SVGs`, y = `Prop Higher Rank`)) +
  geom_point() +
  #add line of best fit
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(aes(label = Dataset), nudge_x = 0.01, nudge_y = 0.01) +
  xlim(0.1, 0.9) +
  theme_bw() +
  labs(x = "Proportion of Unweighted SVGs", y = "Proportion of Small Lengthscale Genes with Higher Rank") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

ggsave(filename=here("plots", "supplementary", "small_lengthscale.png"), 
       plot = fig1, 
       width = 7, 
       height = 7,
       units = "in")
