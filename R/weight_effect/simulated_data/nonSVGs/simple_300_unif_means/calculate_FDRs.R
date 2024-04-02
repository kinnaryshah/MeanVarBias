library(SpatialExperiment)
library(dplyr)
library(ggplot2)

spe_unweighted <- readRDS("spe_nnSVG.rds")
spe_weighted <- readRDS("spe_weighted_nnSVG.rds")

# proportion of rowData(spe_unweighted)$ground_truth_sigma.sq == 0
prop <- sum(rowData(spe_unweighted)$ground_truth_sigma.sq == 0) / dim(spe_unweighted)[1]
# prop = 0.6

# create a loop from the FDR calculation for varying alpha levels from 0.01 to 0.1
# then plot the FDRs for unweighted and weighted models
alpha_levels <- seq(0.01, 0.1, by = 0.01)
fdr_unweighted <- numeric(length(alpha_levels))
fdr_weighted <- numeric(length(alpha_levels))

for (i in 1:length(alpha_levels)) {
  alpha <- alpha_levels[i]
  
  # UNWEIGHTED
  
  rowData(spe_unweighted)$pos_test = rowData(spe_unweighted)$padj < alpha
  rowData(spe_unweighted)$SVG = rowData(spe_unweighted)$ground_truth_sigma.sq != 0
  
  total_genes_unweighted <- dim(spe_unweighted)[1]
  false_positives_unweighted <- sum(rowData(spe_unweighted)$pos_test & !rowData(spe_unweighted)$SVG)
  
  fdr_unweighted[i] <- false_positives_unweighted / sum(rowData(spe_unweighted)$pos_test)
  
  # WEIGHTED
  
  pvalues_weighted <- 1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df = 2)
  pvalues_weighted <- p.adjust(pvalues_weighted, method = "BH")
  rowData(spe_weighted)$pos_test = pvalues_weighted < alpha
  rowData(spe_weighted)$SVG = rowData(spe_weighted)$ground_truth_sigma.sq != 0
  
  total_genes_weighted <- dim(spe_weighted)[1]
  false_positives_weighted <- sum(rowData(spe_weighted)$pos_test & !rowData(spe_weighted)$SVG)
  
  fdr_weighted[i] <- false_positives_weighted / sum(rowData(spe_weighted)$pos_test)
}

# plot using ggplot2
df <- data.frame(alpha = alpha_levels, fdr_unweighted = fdr_unweighted, fdr_weighted = fdr_weighted)
df <- tidyr::gather(df, model, FDR, -alpha)

ggplot(df, aes(x = alpha, y = FDR, color = model)) +
  geom_line() +
  geom_point() +
  labs(title = "FDR vs. Alpha Levels",
       x = "Alpha Level",
       y = "False Discovery Rate",
       color = "Model") +
  theme_bw()

ggsave("FDR_vs_alpha.png")
