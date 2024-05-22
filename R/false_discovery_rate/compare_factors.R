library(ggplot2)
library(ggridges)
library(patchwork)
library(SpatialExperiment)
library(dplyr)
library(ggpubr)

create_FDR_df <- function(spe_unweighted, spe_weighted) {
  
  alpha_levels <- seq(0.001, 0.1, by = 0.001)  
  
  # UNWEIGHTED
  
  
  adj_p_values_unweighted <- rowData(spe_unweighted)$padj
  
  fdr_percent_unweighted <- sapply(alpha_levels, function(alpha) {
    num_positive_tests <- sum(adj_p_values_unweighted <= alpha)
    # calculate the number of false positives as the number of non DE genes with adjusted p-value <= alpha
    num_false_positives <- sum(rowData(spe_unweighted)$padj[rowData(spe_unweighted)$ground_truth_sigma.sq == 0] <= alpha)
    if (num_positive_tests > 0) {
      fdr <- num_false_positives / num_positive_tests 
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
      fdr <- num_false_positives / num_positive_tests
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

create_ridge_plot_weighted <- function(spe_unweighted, spe_weighted) {
  
  #ridge plots separated by noise and signal for weighted with unweighted mean
  frac <- round(dim(spe_weighted)[1]*0.1*0.1)
  n_genes <- dim(spe_weighted)[1]
  
  df <- data.frame(
    rank = rowData(spe_weighted)$weighted_rank,
    mean = rowData(spe_unweighted)$mean
  ) %>% mutate(quantile = findInterval(mean, 
                                       quantile(mean, probs=0:9/10))) %>%
    tibble::rownames_to_column()
  
  # create a new column called "range" with string values
  # this new column contains the range of "mean" values for each "quantile" value
  df <- df %>% 
    group_by(quantile) %>% 
    mutate(range = paste0(round(min(mean),2), " - ", round(max(mean),2))) %>% 
    ungroup()
  
  df_signal <- df %>%
    mutate(quantile = as.factor(quantile)) %>%
    group_by(quantile) %>%
    slice_min(order_by = rank, n = frac) %>%
    mutate(grp = "signal")
  
  indices <- as.integer(df_signal$rowname)
  
  df_background <- df[-indices,] %>%
    mutate(quantile = as.factor(quantile)) %>%
    mutate(grp = "background")
  
  df <- rbind(df_signal, df_background)
  
  # summarize the means of the ranks of "signal" groups within each "quantile"
  # and the means of the ranks of "background" groups within each "quantile"
  # and store the results in a new data frame
  df_rank_means_weighted <- df %>%
    group_by(quantile, grp) %>%
    summarize(mean = mean(rank)) %>%
    ungroup()
  
  ridge_plot <- ggplot(df, aes(x = rank, y = range)) +
    geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3,
                         stat = "binline", bins = 75, draw_baseline = FALSE) +
    theme_ridges(grid = TRUE) +
    labs(
      y = "decile - unw mean of logcounts",
      x = "rank",
      title = ""
    ) +
    coord_cartesian(xlim = c(1, n_genes)) +
    theme_bw() +
    theme(legend.position="none") 
  
  return(ridge_plot)
  
}

create_ridge_plot_unweighted <- function(spe_unweighted) {
  
  frac <- round(dim(spe_unweighted)[1]*0.1*0.1)
  n_genes <- dim(spe_unweighted)[1]
  
  #ridge plots separated by noise and signal for unweighted
  df <- data.frame(
    rank = rowData(spe_unweighted)$rank,
    mean = rowData(spe_unweighted)$mean
  ) %>% mutate(quantile = findInterval(mean, 
                                       quantile(mean, probs=0:9/10))) %>%
    tibble::rownames_to_column()
  
  # create a new column called "range" with string values
  # this new column contains the range of "mean" values for each "quantile" value
  df <- df %>% 
    group_by(quantile) %>% 
    mutate(range = paste0(round(min(mean),2), " - ", round(max(mean),2))) %>% 
    ungroup()
  
  df_signal <- df %>%
    mutate(quantile = as.factor(quantile)) %>%
    group_by(quantile) %>%
    slice_min(order_by = rank, n = frac) %>%
    mutate(grp = "signal")
  
  indices <- as.integer(df_signal$rowname)
  
  df_background <- df[-indices,] %>%
    mutate(quantile = as.factor(quantile)) %>%
    mutate(grp = "background")
  
  df <- rbind(df_signal, df_background)
  
  # summarize the means of the ranks of "signal" groups within each "quantile"
  # and the means of the ranks of "background" groups within each "quantile"
  # and store the results in a new data frame
  df_rank_means_unweighted <- df %>%
    group_by(quantile, grp) %>%
    summarize(mean = mean(rank)) %>%
    ungroup()
  
  ridge_plot <- ggplot(df, aes(x = rank, y = range)) +
    geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3,
                         stat = "binline", bins = 75, draw_baseline = FALSE) +
    theme_ridges(grid = TRUE) +
    labs(
      y = "decile - unw mean of logcounts",
      x = "rank",
      title = ""
    ) +
    theme_bw() + 
    theme(legend.position="none") +
    coord_cartesian(xlim = c(1, n_genes)) 

  return(ridge_plot)
}


# big picture
# 4 rows of plots
# first row: ridge plot for unweighted
# second row: ridge plot for weighted
# third row: FDR for unweighted and weighted
# fourth row: TPR and TNR for unweighted and weighted

# 5 columns of plots
# first column: best
# second column: high_genes
# third column: high_sigma
# fourth column: high_beta
# fifth column: high_spots

# create a function to make a list of plots for each column
# inputs should be file directory
# outputs should be a list of plots

create_col_plots <- function(file_dir) {
  spe_unweighted <- readRDS(paste0(file_dir, "/spe_nnSVG_1.rds"))
  spe_weighted <- readRDS(paste0(file_dir, "/spe_weighted_nnSVG_1.rds"))
  
  # ridge plot for unweighted
  p_ridge_unweighted <- create_ridge_plot_unweighted(spe_unweighted)
  
  # ridge plot for weighted
  p_ridge_weighted <- create_ridge_plot_weighted(spe_unweighted, spe_weighted)
  
  # FDR for unweighted and weighted
  df <- data.frame()
  for (i in 1:10) {
    print(i)
    spe_unweighted <- readRDS(paste0(file_dir, "/spe_nnSVG_", i, ".rds"))
    spe_weighted <- readRDS(paste0(file_dir, "/spe_weighted_nnSVG_", i, ".rds"))
    
    df_temp <- create_FDR_df(spe_unweighted, spe_weighted)
    df_temp$sim <- i
    
    df <- rbind(df, df_temp)
  }
  
  df_avg <- df %>%
    group_by(alpha) %>%
    summarise(fdr_percent_unweighted_sd = sd(fdr_percent_unweighted),
              fdr_percent_weighted_sd = sd(fdr_percent_weighted),
              fdr_percent_unweighted = mean(fdr_percent_unweighted),
              fdr_percent_weighted = mean(fdr_percent_weighted))
  
  p_FDR <- df_avg %>%
    ggplot(aes(x = alpha)) +
    geom_line(aes(y = fdr_percent_unweighted), color = "blue") +
    geom_line(aes(y = fdr_percent_weighted), color = "red") +
    geom_ribbon(aes(ymin = fdr_percent_unweighted - fdr_percent_unweighted_sd, ymax = fdr_percent_unweighted + fdr_percent_unweighted_sd), fill = "blue") +
    geom_ribbon(aes(ymin = fdr_percent_weighted - fdr_percent_weighted_sd, ymax = fdr_percent_weighted + fdr_percent_weighted_sd), fill = "red") +
    ylim(0, 0.4) +
    labs(title = "", x = "alpha", y = "FDR") +
    theme(legend.position="none") +
    theme_bw() 
  
  # TPR and TNR for unweighted and weighted
  df <- data.frame()
  for (i in 1:10) {
    print(i)
    spe_unweighted <- readRDS(paste0(file_dir, "/spe_nnSVG_", i, ".rds"))
    spe_weighted <- readRDS(paste0(file_dir, "/spe_weighted_nnSVG_", i, ".rds"))
    
    df_temp <- create_TNR_TPR_df(spe_unweighted, spe_weighted)
    df_temp$sim <- i
    
    df <- rbind(df, df_temp)
  }
  
  df_avg <- df %>%
    group_by(fdr_levels) %>%
    summarise(tpr_levels_unweighted_sd = sd(tpr_levels_unweighted),
              tnr_levels_unweighted_sd = sd(tnr_levels_unweighted),
              tpr_levels_weighted_sd = sd(tpr_levels_weighted),
              tnr_levels_weighted_sd = sd(tnr_levels_weighted),
              tpr_levels_unweighted = mean(tpr_levels_unweighted),
              tnr_levels_unweighted = mean(tnr_levels_unweighted),
              tpr_levels_weighted = mean(tpr_levels_weighted),
              tnr_levels_weighted = mean(tnr_levels_weighted))
  
  p_TNR <- df_avg %>%
    ggplot(aes(x = fdr_levels)) +
    geom_line(aes(y = tnr_levels_unweighted), color = "blue", linetype = "dashed") +
    geom_line(aes(y = tnr_levels_weighted), color = "red", linetype = "dashed") +
    geom_ribbon(aes(ymin = tnr_levels_unweighted - tnr_levels_unweighted_sd, ymax = tnr_levels_unweighted + tnr_levels_unweighted_sd), fill = "blue") +
    geom_ribbon(aes(ymin = tnr_levels_weighted - tnr_levels_weighted_sd, ymax = tnr_levels_weighted + tnr_levels_weighted_sd), fill = "red") +
    labs(title = "", x = "FDR", y = "TNR") +
    scale_linetype_manual(values = c("dashed", "solid")) +
    ylim(0, 1) +
    theme_bw()
  
  p_TPR <- df_avg %>%
    ggplot(aes(x = fdr_levels)) +
    geom_line(aes(y = tpr_levels_unweighted), color = "blue") +
    geom_line(aes(y = tpr_levels_weighted), color = "red") +
    geom_ribbon(aes(ymin = tpr_levels_unweighted - tpr_levels_unweighted_sd, ymax = tpr_levels_unweighted + tpr_levels_unweighted_sd), fill = "blue") +
    geom_ribbon(aes(ymin = tpr_levels_weighted - tpr_levels_weighted_sd, ymax = tpr_levels_weighted + tpr_levels_weighted_sd), fill = "red") +
    labs(title = "", x = "FDR", y = "TPR") +
    scale_linetype_manual(values = c("dashed", "solid")) +
    ylim(0, 1) +
    theme_bw()
  
  return(list(p_ridge_unweighted, p_ridge_weighted, p_FDR, p_TNR, p_TPR))
}


percent <- "70"

best <- "reps_high_range_sigma.sq_high_beta_1000_genes_SVG_percentages_subset"
high_genes <- "reps_high_range_sigma.sq_high_beta_3000_genes_SVG_percentages_subset"
high_sigma <- "reps_high_sigma.sq_high_beta_1000_genes_SVG_percentages_subset"
high_beta <- "reps_high_range_sigma.sq_huge_beta_1000_genes_SVG_percentages_subset"
high_spots <- "reps_high_range_sigma.sq_high_beta_1000_genes_2000_spots_SVG_percentages_subset"

best_directory <- paste0(best,"/",percent,"_percent_nonSVGs")
best_list <- create_col_plots(best_directory)

high_genes_directory <- paste0(high_genes,"/",percent,"_percent_nonSVGs")
high_genes_list <- create_col_plots(high_genes_directory)

high_spots_directory <- paste0(high_spots,"/",percent,"_percent_nonSVGs")
high_spots_list <- create_col_plots(high_spots_directory)

#plot <- ggarrange(plotlist = c(best_list, high_genes_list, high_sigma_list, high_beta_list, high_spots_list),
#                  ncol = 4, nrow = 5)

#pdf(paste0(percent,"_percent_plots.pdf"), width = 10, height = 15)
#annotate_figure(plot, top = text_grob(paste0(percent, "% Non-SVGs"), face = "bold", size = 14),
#                left = "Increase Spots                                          High Beta                                                    High Sigma                                          Increase Genes                                                             Best")
#dev.off()

plot <- ggarrange(plotlist = c(best_list, high_genes_list, high_spots_list),
                  ncol = 5, nrow = 3)

pdf(paste0(percent,"_percent_plots.pdf"), width = 10, height = 15)
annotate_figure(plot, top = text_grob(paste0(percent, "% Non-SVGs"), face = "bold", size = 14),
left = "Increase Spots                                                                                               Increase Genes                                                                                                 Best")
dev.off()

