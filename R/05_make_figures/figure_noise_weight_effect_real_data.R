library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(viridis)
library(ggpubr)
library(purrr)
library(patchwork)
library(GGally)
library(ggridges)

# overall idea: create a function to make four plots for each dataset
# then, use the function for each dataset of interest and combine the plots

ridge_plots <- function(spe_unweighted, spe_weighted, dataset_name){
  
  n_genes <- dim(spe_weighted)[1]
  frac <- round(dim(spe_weighted)[1]*0.1*0.1)
  
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
  
  unweighted_ridge_combined <- df %>% 
    mutate(quantile = as.factor(quantile)) |> 
    ggplot() +
    geom_density_ridges2(aes(x = rank, y = range), rel_min_height = 0.02) +
    theme_ridges(grid = TRUE) +
    labs(
      y = paste0(dataset_name, " quantile - unw mean of logcounts"),
      x = "rank",
      title = "unweighted"
    ) +
    coord_cartesian(xlim = c(1, n_genes)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
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
  
  unweighted_ridge <- ggplot(df, aes(x = rank, y = range)) +
    geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3,
                         stat = "binline", bins = 75, draw_baseline = FALSE) +
    theme_ridges(grid = TRUE) +
    labs(
      y = "",
      x = "rank",
      title = ""
    ) +
    guides(fill=guide_legend(title="group")) +
    coord_cartesian(xlim = c(1, n_genes)) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  #ridge plots separated by noise and signal for weighted with unweighted mean
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
  
  weighted_ridge_combined_unw_mean <- df %>% 
    mutate(quantile = as.factor(quantile)) |> 
    ggplot() +
    geom_density_ridges2(aes(x = rank, y = range), rel_min_height = 0.02) +
    theme_ridges(grid = TRUE) +
    labs(
      y = "",
      x = "rank",
      title = "weighted"
    ) +
    coord_cartesian(xlim = c(1, n_genes)) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
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
  
  weighted_ridge_unw_mean <- ggplot(df, aes(x = rank, y = range)) +
    geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3,
                         stat = "binline", bins = 75, draw_baseline = FALSE) +
    theme_ridges(grid = TRUE) +
    labs(
      y = "",
      x = "rank",
      title = ""
    ) +
    guides(fill=guide_legend(title="group")) +
    coord_cartesian(xlim = c(1, n_genes)) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  return(list(unweighted_ridge_combined, unweighted_ridge, weighted_ridge_combined_unw_mean, weighted_ridge_unw_mean))
}

spe_unweighted <- readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG.rds"))
spe_weighted <- readRDS(here("outputs", "results", "spe_humanDLPFC_noise_weighted_nnSVG.rds"))
plots_DLPFC <- ridge_plots(spe_unweighted, spe_weighted, "DLPFC")

spe_unweighted <- readRDS(here("outputs", "results", "spe_humanOvarian_nnSVG.rds"))
spe_weighted <- readRDS(here("outputs", "results", "spe_humanOvarian_noise_weighted_nnSVG.rds"))
plots_ovarian <- ridge_plots(spe_unweighted, spe_weighted, "Ovarian")

spe_unweighted <- readRDS(here("outputs", "results", "spe_humanLC_nnSVG.rds"))
spe_weighted <- readRDS(here("outputs", "results", "spe_humanLC_noise_weighted_nnSVG.rds"))
plots_LC <- ridge_plots(spe_unweighted, spe_weighted, "LC")

spe_unweighted <- readRDS(here("outputs", "results", "spe_humanHPC_nnSVG.rds"))
spe_weighted <- readRDS(here("outputs", "results", "spe_humanHPC_noise_weighted_nnSVG.rds"))
plots_HPC <- ridge_plots(spe_unweighted, spe_weighted, "HPC")

spe_unweighted <- readRDS(here("outputs", "results", "spe_humanBreast_nnSVG.rds"))
spe_weighted <- readRDS(here("outputs", "results", "spe_humanBreast_noise_weighted_nnSVG.rds"))
plots_breast <- ridge_plots(spe_unweighted, spe_weighted, "Breast")

spe_unweighted <- readRDS(here("outputs", "results", "spe_humanLobularBreast_nnSVG.rds"))
spe_weighted <- readRDS(here("outputs", "results", "spe_humanLobularBreast_noise_weighted_nnSVG.rds"))
plots_lobular_breast <- ridge_plots(spe_unweighted, spe_weighted, "Lobular Breast")

pdf(here("plots", "main", "noise_weight_effect_real_data.pdf"))
wrap_plots(c(plots_DLPFC, plots_ovarian), guides="collect",
           ncol=4, nrow=2)
wrap_plots(c(plots_LC, plots_HPC), guides="collect",
           ncol=4, nrow=2)
wrap_plots(c(plots_breast, plots_lobular_breast), guides="collect",
           ncol=4, nrow=2)
dev.off()
