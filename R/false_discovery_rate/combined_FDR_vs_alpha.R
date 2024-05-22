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
  print(i)
  spe_unweighted <- readRDS(paste0("spe_nnSVG_", i, ".rds"))
  
  p1 <- data.frame(pval = rowData(spe_unweighted)$pval) %>%
    ggplot(aes(x = pval)) +
    geom_histogram(binwidth = 0.01) +
    ggtitle(paste0("unweighted p-value distribution for ", i, "0% nonSVGs")) +
    # add limits to y-axis
    ylim(0, 1000) +
    theme_minimal()
  
  spe_weighted <- readRDS(paste0("spe_weighted_nnSVG_", i, ".rds"))
  p2 <- data.frame(pval = 1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df=2)) %>%
    ggplot(aes(x = pval)) +
    geom_histogram(binwidth = 0.01) +
    ggtitle(paste0("weighted p-value distribution for ", i, "0% nonSVGs")) +
    # add limits to y-axis
    ylim(0, 1000) +
    theme_minimal()
  
  # save plots
  pdf(paste0("p_value_distribution_", i, ".pdf"))
  print(p1)
  print(p2)
  dev.off()
  
}

# compute Kolmogorov–Smirnov (K-S) distance between the distribution of the computed p-values and the uniform distribution
# for each replicate
ks_unweighted <- c()
ks_weighted <- c()

for (i in 3:7) {
  spe_unweighted <- readRDS(paste0("spe_nnSVG_", i, ".rds"))
  ks_unweighted <- c(ks_unweighted, ks.test(rowData(spe_unweighted)$pval, "punif")$statistic)
  
  spe_weighted <- readRDS(paste0("spe_weighted_nnSVG_", i, ".rds"))
  ks_weighted <- c(ks_weighted, ks.test(1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df=2), "punif")$statistic)
}

# 2. plot average FDR for different alpha values
# one plot combining all replicates
# make a function with input of spe_unweighted and spe_weighted, return a data frame of fdr_unweighted and fdr_weighted 
create_FDR_df <- function(spe_unweighted, spe_weighted) {
  
  alpha_levels <- seq(0.001, 0.1, by = 0.001)  
  
  # UNWEIGHTED
  
  
  adj_p_values_unweighted <- rowData(spe_unweighted)$padj
  
  fdr_percent_unweighted <- sapply(alpha_levels, function(alpha) {
    num_positive_tests <- sum(adj_p_values_unweighted <= alpha)
    # calculate the number of false positives as the number of non DE genes with adjusted p-value <= alpha
    num_false_positives <- sum(rowData(spe_unweighted)$padj[rowData(spe_unweighted)$ground_truth_sigma.sq == 0] <= alpha)
    if (num_positive_tests > 0) {
      fdr <- num_false_positives / num_positive_tests * 100
    } else {
      fdr <- 0
    }
    return(fdr)
  })
  
  
  # WEIGHTED
  
  adj_p_values_weighted <- p.adjust(1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df=2), method = "BH")
  rowData(spe_weighted)$weighted_padj <- adj_p_values_weighted
  
  fdr_percent_weighted <- sapply(alpha_levels, function(alpha) {
    num_positive_tests <- sum(adj_p_values_weighted <= alpha, na.rm = TRUE)
    # calculate the number of false positives as the number of non DE genes with adjusted p-value <= alpha
    num_false_positives <- sum(rowData(spe_weighted)$weighted_padj[rowData(spe_weighted)$ground_truth_sigma.sq == 0] <= alpha, na.rm = TRUE)
    if (num_positive_tests > 0) {
      fdr <- num_false_positives / num_positive_tests * 100
    } else {
      fdr <- 0
    }
    return(fdr)
  })
  
  df <- data.frame(fdr_percent_unweighted = fdr_percent_unweighted,
                   fdr_percent_weighted = fdr_percent_weighted,
                   alpha = alpha_levels)
  
  return(df)
}


# use the function to create the data frame for the files from all the folders from fdr_1 to fdr_8
# create a column for the fdr_1 to fdr_8 

df <- data.frame()
for (i in 1:5) {
  print(i)
  spe_unweighted <- readRDS(paste0("spe_nnSVG_", i, ".rds"))
  spe_weighted <- readRDS(paste0("spe_weighted_nnSVG_", i, ".rds"))
  
  df_temp <- create_FDR_df(spe_unweighted, spe_weighted)
  df_temp$sim <- i
  
  df <- rbind(df, df_temp)
}

# average FDR for different alpha values
df_avg <- df %>%
  group_by(alpha) %>%
  summarise(fdr_percent_unweighted = mean(fdr_percent_unweighted),
            fdr_percent_weighted = mean(fdr_percent_weighted))

# plot average FDR for different alpha values
p_FDR <- df_avg %>%
  ggplot(aes(x = alpha)) +
  geom_line(aes(y = fdr_percent_unweighted, color = "unweighted")) +
  geom_line(aes(y = fdr_percent_weighted, color = "weighted")) +
  ggtitle("Average FDR vs alpha") +
  theme_bw() +
  theme(legend.position = "top")

# 3. plot average TPR and TPR vs alpha
# one plot combining all replicates
create_TNR_TPR_df <- function(spe_unweighted, spe_weighted) {
  
  fdr_levels <- seq(0.001, 0.1, by = 0.001)  
  tpr_levels_unweighted <- numeric(length(fdr_levels))
  tnr_levels_unweighted <- numeric(length(fdr_levels))
  tpr_levels_weighted <- numeric(length(fdr_levels))
  tnr_levels_weighted <- numeric(length(fdr_levels))
  
  # UNWEIGHTED
  
  for (k in 1:length(fdr_levels)) {
    tp <- sum(rowData(spe_unweighted)$ground_truth_sigma.sq != 0 & rowData(spe_unweighted)$padj <= fdr_levels[k])
    tn <- sum(rowData(spe_unweighted)$ground_truth_sigma.sq == 0 & rowData(spe_unweighted)$padj > fdr_levels[k])
    fn <- sum(rowData(spe_unweighted)$ground_truth_sigma.sq != 0 & rowData(spe_unweighted)$padj > fdr_levels[k])
    fp <- sum(rowData(spe_unweighted)$ground_truth_sigma.sq == 0 & rowData(spe_unweighted)$padj <= fdr_levels[k])
    
    tpr_levels_unweighted[k] <- tp / (tp + fn)
    tnr_levels_unweighted[k] <- tn / (tn + fp)
  }
  
  # WEIGHTED
  
  adj_p_values_weighted <- p.adjust(1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df=2), method = "BH")
  rowData(spe_weighted)$weighted_padj <- adj_p_values_weighted
  
  for (k in 1:length(fdr_levels)) {
    tp <- sum(rowData(spe_weighted)$ground_truth_sigma.sq != 0 & rowData(spe_weighted)$weighted_padj <= fdr_levels[k], na.rm = TRUE)
    tn <- sum(rowData(spe_weighted)$ground_truth_sigma.sq == 0 & rowData(spe_weighted)$weighted_padj > fdr_levels[k], na.rm = TRUE)
    fn <- sum(rowData(spe_weighted)$ground_truth_sigma.sq != 0 & rowData(spe_weighted)$weighted_padj > fdr_levels[k], na.rm = TRUE)
    fp <- sum(rowData(spe_weighted)$ground_truth_sigma.sq == 0 & rowData(spe_weighted)$weighted_padj <= fdr_levels[k], na.rm = TRUE)
    
    tpr_levels_weighted[k] <- tp / (tp + fn)
    tnr_levels_weighted[k] <- tn / (tn + fp)
  }
  
  df <- data.frame(tpr_levels_unweighted = tpr_levels_unweighted,
                   tnr_levels_unweighted = tnr_levels_unweighted,
                   tpr_levels_weighted = tpr_levels_weighted,
                   tnr_levels_weighted = tnr_levels_weighted,
                   fdr_levels = fdr_levels)
  
  return(df)
}

df <- data.frame()
for (i in 1:5) {
  print(i)
  spe_unweighted <- readRDS(paste0("spe_nnSVG_", i, ".rds"))
  spe_weighted <- readRDS(paste0("spe_weighted_nnSVG_", i, ".rds"))
  
  df_temp <- create_TNR_TPR_df(spe_unweighted, spe_weighted)
  df_temp$sim <- i
  
  df <- rbind(df, df_temp)
}

# average TPR and TNR for different FDR levels
df_avg <- df %>%
  group_by(fdr_levels) %>%
  summarise(tpr_levels_unweighted = mean(tpr_levels_unweighted),
            tnr_levels_unweighted = mean(tnr_levels_unweighted),
            tpr_levels_weighted = mean(tpr_levels_weighted),
            tnr_levels_weighted = mean(tnr_levels_weighted))

# plot average TPR and TNR for different FDR levels
p_TPR <- df_avg %>%
  ggplot(aes(x = fdr_levels)) +
  geom_line(aes(y = tpr_levels_unweighted, color = "unweighted")) +
  geom_line(aes(y = tpr_levels_weighted, color = "weighted")) +
  ggtitle("Average TPR vs FDR") +
  theme_bw() +
  theme(legend.position = "top")

p_TNR <- df_avg %>%
  ggplot(aes(x = fdr_levels)) +
  geom_line(aes(y = tnr_levels_unweighted, color = "unweighted")) +
  geom_line(aes(y = tnr_levels_weighted, color = "weighted")) +
  ggtitle("Average TNR vs FDR") +
  theme_bw() +
  theme(legend.position = "top")

# save plots
pdf("FDR_TPR_TNR.pdf")
print(p_FDR)
print(p_TPR)
print(p_TNR)
dev.off()
