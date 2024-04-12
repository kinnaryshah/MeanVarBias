library(SpatialExperiment)
library(dplyr)
library(ggplot2)

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
for (i in 1:8) {
  print(i)
  spe_unweighted <- readRDS(paste0("fdr_", i, "/spe_nnSVG.rds"))
  spe_weighted <- readRDS(paste0("fdr_", i, "/spe_weighted_nnSVG.rds"))
  
  df_temp <- create_FDR_df(spe_unweighted, spe_weighted)
  df_temp$sim <- i+1
  
  df <- rbind(df, df_temp)
}

# add another df_temp to df using i=1 from the folder "../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/"
spe_unweighted <- readRDS("../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/spe_nnSVG.rds")
spe_weighted <- readRDS("../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/spe_weighted_nnSVG.rds")

df_temp <- create_FDR_df(spe_unweighted, spe_weighted)
df_temp$sim <- 1

df <- rbind(df, df_temp)

# plot the FDRs for unweighted and weighted models
# separate plots for each fdr 1 to 9
# use faceting from ggplot2
# x axis is the alpha level
# y axis is the FDR percentage
# color the lines by the model type (unweighted or weighted)

p <- ggplot(df, aes(x = alpha, y = fdr_percent_unweighted, color = "unweighted")) +
  geom_line() +
  geom_line(aes(y = fdr_percent_weighted, color = "weighted")) +
  facet_wrap(~sim) +
  labs(title = "FDR vs significance level",
       x = "Significance level",
       y = "FDR (%)",
       color = "Model") +
  theme_bw()

# save the plot
ggsave("combined_FDR_vs_alpha.png", p, width = 10, height = 10)

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
for (i in 1:8) {
  print(i)
  spe_unweighted <- readRDS(paste0("fdr_", i, "/spe_nnSVG.rds"))
  spe_weighted <- readRDS(paste0("fdr_", i, "/spe_weighted_nnSVG.rds"))
  
  df_temp <- create_TNR_TPR_df(spe_unweighted, spe_weighted)
  df_temp$sim <- i+1
  
  df <- rbind(df, df_temp)
}

# add another df_temp to df using i=1 from the folder "../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/"
spe_unweighted <- readRDS("../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/spe_nnSVG.rds")
spe_weighted <- readRDS("../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/spe_weighted_nnSVG.rds")

df_temp <- create_TNR_TPR_df(spe_unweighted, spe_weighted)
df_temp$sim <- 1

df <- rbind(df, df_temp)

# each plot should have two lines (one for unweighted and one for weighted)
# x axis is the FDR level
# y axis is the TNR or TPR
# color the lines by the model type (unweighted or weighted)
# dash or solid lines for TNR and TPR

p <- ggplot(df, aes(x = fdr_levels)) +
  geom_line(aes(y = tnr_levels_unweighted, color = "unweighted", linetype = "dashed")) +
  geom_line(aes(y = tnr_levels_weighted, color = "weighted", linetype = "dashed")) +
  facet_wrap(~sim) +
  labs(title = "TNR vs FDR",
       x = "FDR",
       y = "TNR",
       color = "Model") +
  theme_bw()

# save the plot
ggsave("combined_TNR.png", p, width = 10, height = 10)

p <- ggplot(df, aes(x = fdr_levels)) +
  geom_line(aes(y = tpr_levels_unweighted, color = "unweighted", linetype = "solid")) +
  geom_line(aes(y = tpr_levels_weighted, color = "weighted", linetype = "solid")) +
  facet_wrap(~sim) +
  labs(title = "TPR vs FDR",
       x = "FDR",
       y = "TPR",
       color = "Model") +
  theme_bw()

# save the plot
ggsave("combined_TPR.png", p, width = 10, height = 10)


create_p_val_df <- function(spe_unweighted, spe_weighted) {
  
  # UNWEIGHTED
  
  # Sort unadjusted p-values and calculate FDR
  sorted_p_values_unweighted <- sort(rowData(spe_unweighted)$pval)
  
  # WEIGHTED
  
  # Sort unadjusted p-values and calculate FDR
  pvalues_weighted <- 1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df = 2)
  sorted_p_values_weighted <- sort(pvalues_weighted, na.last = TRUE)
  
  df <- data.frame(sorted_p_values_weighted = sorted_p_values_weighted,
                   sorted_p_values_unweighted = sorted_p_values_unweighted)
  
  return(df)
}


# use the function to create the data frame for the files from all the folders from fdr_1 to fdr_8
# create a column for the fdr_1 to fdr_8 

df <- data.frame()
for (i in 1:8) {
  print(i)
  spe_unweighted <- readRDS(paste0("fdr_", i, "/spe_nnSVG.rds"))
  spe_weighted <- readRDS(paste0("fdr_", i, "/spe_weighted_nnSVG.rds"))
  
  df_temp <- create_p_val_df(spe_unweighted, spe_weighted)
  df_temp$sim <- i+1
  
  df <- rbind(df, df_temp)
}

# add another df_temp to df using i=1 from the folder "../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/"
spe_unweighted <- readRDS("../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/spe_nnSVG.rds")
spe_weighted <- readRDS("../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/spe_weighted_nnSVG.rds")

df_temp <- create_p_val_df(spe_unweighted, spe_weighted)
df_temp$sim <- 1

df <- rbind(df, df_temp)

# plot p value distributions for sim 1 only
# two separate plots for weighted and unweighted
# x axis is the p value
# y axis is the density
# color the lines by the model type (unweighted or weighted)
sub_df <- df[df$sim == 1, ]

p1 <- ggplot(sub_df, aes(x = sorted_p_values_unweighted)) +
  geom_histogram(aes(y = ..density..), bins = 50, color = "black", fill = "blue", alpha = 0.5) +
  labs(title = "P-value distribution for unweighted model",
       x = "P-value",
       y = "Density") +
  theme_bw()

p2 <- ggplot(sub_df, aes(x = sorted_p_values_weighted)) +
  geom_histogram(aes(y = ..density..), bins = 50, color = "black", fill = "blue", alpha = 0.5) +
  labs(title = "P-value distribution for weighted model",
       x = "P-value",
       y = "Density") +
  theme_bw()

# save the plots together
pdf("p_val_distributions_sim1.pdf")
print(p1)
print(p2)
dev.off()

# plot p value distributions for sim 9 only
# two separate plots for weighted and unweighted
# x axis is the p value
# y axis is the density
# color the lines by the model type (unweighted or weighted)
sub_df <- df[df$sim == 9, ]

p1 <- ggplot(sub_df, aes(x = sorted_p_values_unweighted)) +
  geom_histogram(aes(y = ..density..), bins = 50, color = "black", fill = "blue", alpha = 0.5) +
  labs(title = "P-value distribution for unweighted model",
       x = "P-value",
       y = "Density") +
  theme_bw()

p2 <- ggplot(sub_df, aes(x = sorted_p_values_weighted)) +
  geom_histogram(aes(y = ..density..), bins = 50, color = "black", fill = "blue", alpha = 0.5) +
  labs(title = "P-value distribution for weighted model",
       x = "P-value",
       y = "Density") +
  theme_bw()

# save the plots together
pdf("p_val_distributions_sim9.pdf")
print(p1)
print(p2)
dev.off()

create_filtered_p_val_df <- function(spe_unweighted, spe_weighted) {
  
  # UNWEIGHTED
  
  # remove lowly expressed genes
  indices <- rowData(spe_unweighted)$mean > 0.75
  spe_unweighted <- spe_unweighted[indices, ]
  sorted_p_values_unweighted <- sort(rowData(spe_unweighted)$pval)
  
  # WEIGHTED
  # remove lowly expressed genes
  # use unweighted mean to filter out lowly expressed genes since this is the true mean
  spe_weighted <- spe_weighted[indices, ]
  pvalues_weighted <- 1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df = 2)
  sorted_p_values_weighted <- sort(pvalues_weighted, na.last = TRUE)
  
  df <- data.frame(sorted_p_values_weighted = sorted_p_values_weighted,
                   sorted_p_values_unweighted = sorted_p_values_unweighted)
  
  return(df)
}


# use the function to create the data frame for the files from all the folders from fdr_1 to fdr_8
# create a column for the fdr_1 to fdr_8 

df <- data.frame()
for (i in 1:8) {
  print(i)
  spe_unweighted <- readRDS(paste0("fdr_", i, "/spe_nnSVG.rds"))
  spe_weighted <- readRDS(paste0("fdr_", i, "/spe_weighted_nnSVG.rds"))
  
  df_temp <- create_filtered_p_val_df(spe_unweighted, spe_weighted)
  df_temp$sim <- i+1
  
  df <- rbind(df, df_temp)
}

# add another df_temp to df using i=1 from the folder "../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/"
spe_unweighted <- readRDS("../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/spe_nnSVG.rds")
spe_weighted <- readRDS("../weight_effect/simulated_data/nonSVGs/simple_300_unif_means/spe_weighted_nnSVG.rds")

df_temp <- create_filtered_p_val_df(spe_unweighted, spe_weighted)
df_temp$sim <- 1

df <- rbind(df, df_temp)

# plot p value distributions for sim 1 only
# two separate plots for weighted and unweighted
# x axis is the p value
# y axis is the density
# color the lines by the model type (unweighted or weighted)
sub_df <- df[df$sim == 1, ]

p1 <- ggplot(sub_df, aes(x = sorted_p_values_unweighted)) +
  geom_histogram(aes(y = ..density..), bins = 50, color = "black", fill = "blue", alpha = 0.5) +
  labs(title = "P-value distribution for unweighted model",
       x = "P-value",
       y = "Density") +
  theme_bw()

p2 <- ggplot(sub_df, aes(x = sorted_p_values_weighted)) +
  geom_histogram(aes(y = ..density..), bins = 50, color = "black", fill = "blue", alpha = 0.5) +
  labs(title = "P-value distribution for weighted model",
       x = "P-value",
       y = "Density") +
  theme_bw()

# save the plots together
pdf("p_val_distributions_sim1_filtered.pdf")
print(p1)
print(p2)
dev.off()

# plot p value distributions for sim 9 only
# two separate plots for weighted and unweighted
# x axis is the p value
# y axis is the density
# color the lines by the model type (unweighted or weighted)
sub_df <- df[df$sim == 9, ]

p1 <- ggplot(sub_df, aes(x = sorted_p_values_unweighted)) +
  geom_histogram(aes(y = ..density..), bins = 50, color = "black", fill = "blue", alpha = 0.5) +
  labs(title = "P-value distribution for unweighted model",
       x = "P-value",
       y = "Density") +
  theme_bw()

p2 <- ggplot(sub_df, aes(x = sorted_p_values_weighted)) +
  geom_histogram(aes(y = ..density..), bins = 50, color = "black", fill = "blue", alpha = 0.5) +
  labs(title = "P-value distribution for weighted model",
       x = "P-value",
       y = "Density") +
  theme_bw()

# save the plots together
pdf("p_val_distributions_sim9_filtered.pdf")
print(p1)
print(p2)
dev.off()



