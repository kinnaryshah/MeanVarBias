library(SpatialExperiment)
library(here)
library(dplyr)
library(ggplot2)


small_lengthscale_effect <- function(spe_unweighted, spe_weighted, low, high) {
  idx <- which(rowData(spe_unweighted)$phi > low & rowData(spe_unweighted)$phi < high)
  
  ranks <- data.frame(
    gene = rowData(spe_unweighted)$gene_name[idx],
    unweighted = rowData(spe_unweighted)$rank[idx],
    weighted = rowData(spe_weighted)$weighted_rank[idx],
    diff = rowData(spe_unweighted)$rank[idx] - rowData(spe_weighted)$weighted_rank[idx]
  )
  
  higher_rank_genes <- ranks$gene[ranks$diff > 0]
  
  prop_unweighted_SVGs <- sum(rowData(spe_unweighted)$padj < 0.05) / dim(spe_unweighted)[1]
  prop_weighted_SVGs <- sum(rowData(spe_weighted)$weighted_padj < 0.05, na.rm = TRUE) / dim(spe_weighted)[1]
  
  print(paste0("Prop of unweighted SVGs: ", round(prop_unweighted_SVGs, 4),
               ", Prop of weighted SVGs: ", round(prop_weighted_SVGs, 4),
               ", Total genes: ", dim(spe_unweighted)[1]))
  
  print("Genes with higher ranks after weighting:")
  print(higher_rank_genes)
  
  return(list(
    higher_rank_count = length(higher_rank_genes),
    total_small_lengthscale = length(idx),
    prop_unweighted_SVGs = prop_unweighted_SVGs,
    prop_weighted_SVGs = prop_weighted_SVGs,
    total_genes = dim(spe_unweighted)[1]
  ))
}

dat <- c("Ovarian", "Breast", "LobularBreast", "SubtypeBreast")

results <- data.frame()
for (d in dat) {
  spe_unweighted <- readRDS(here("outputs", "results", paste0("spe_human", d, "_nnSVG.rds")))
  spe_weighted <- readRDS(here("outputs", "results", paste0("spe_human", d, "_weighted_nnSVG.rds")))
  
  res <- small_lengthscale_effect(spe_unweighted, spe_weighted, 0.011, 0.025)
  results <- rbind(results, c(d, res$higher_rank_count, res$total_small_lengthscale,
                              res$prop_unweighted_SVGs, res$prop_weighted_SVGs, res$total_genes))
}

colnames(results) <- c("Dataset", "Higher Rank", "Total Small Lengthscale", 
                       "Prop Unweighted SVGs", "Prop Weighted SVGs", "Total Genes")
results$`Prop Higher Rank` <- results$`Higher Rank` / results$`Total Small Lengthscale`

results <- results[, c(1, 2, 3, 7, 4, 5, 6)]

print(results)