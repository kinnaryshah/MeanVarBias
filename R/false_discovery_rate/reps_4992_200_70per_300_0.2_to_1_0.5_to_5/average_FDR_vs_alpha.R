library(SpatialExperiment)
library(dplyr)
library(ggplot2)

# 1. visualize p-value distributions for 5 replicates
# 2. plot average FDR for different alpha values
# 3. plot average TPR vs alpha
# 4. plot average TNR vs alpha

# 1. visualize p-value distributions for 5 replicates
# separate plots 

for (i in 1:5) {
  spe_unweighted <- readRDS(paste0("spe_nnSVG_", i, ".rds"))
  
  p1 <- data.frame(pval = rowData(spe_unweighted)$pval) %>%
    ggplot(aes(x = pval)) +
    geom_histogram(binwidth = 0.01) +
    ggtitle(paste0("unweighted p-value distribution for replicate ", i)) +
    theme_minimal()

  spe_weighted <- readRDS(paste0("spe_weighted_nnSVG_", i, ".rds"))
  p2 <- data.frame(pval = 1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df=2)) %>%
    ggplot(aes(x = pval)) +
    geom_histogram(binwidth = 0.01) +
    ggtitle(paste0("weighted p-value distribution for replicate ", i)) +
    theme_minimal()
  
  # save plots
  pdf(paste0("p_value_distribution_", i, ".pdf"))
  print(p1)
  print(p2)
  dev.off()
    
}

for (i in 1:5) {
  spe_unweighted <- readRDS(paste0("spe_nnSVG_", i, ".rds"))
  indices <- rowData(spe_unweighted)$mean > 0.75
  spe_unweighted <- spe_unweighted[indices, ]
  
  p1 <- data.frame(pval = rowData(spe_unweighted)$pval) %>%
    ggplot(aes(x = pval)) +
    geom_histogram(binwidth = 0.01) +
    ggtitle(paste0("unweighted p-value distribution for replicate ", i)) +
    theme_minimal()
 
  spe_weighted <- readRDS(paste0("spe_weighted_nnSVG_", i, ".rds"))
  spe_weighted <- spe_weighted[indices, ]
  p2 <- data.frame(pval = 1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df=2)) %>%
    ggplot(aes(x = pval)) +
    geom_histogram(binwidth = 0.01) +
    ggtitle(paste0("weighted p-value distribution for replicate ", i)) +
    theme_minimal()
  
  # save plots
  pdf(paste0("p_value_distribution_filtered_", i, ".pdf"))
  print(p1)
  print(p2)
  dev.off()
  
}

