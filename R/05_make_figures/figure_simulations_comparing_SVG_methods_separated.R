library(ggplot2)
library(SpatialExperiment)
library(patchwork)
library(GGally)
library(dplyr)
library(ggridges)
library(here)

sim_dir <- "reps_968_50_50per_1000_0.2_to_3_0.5_to_9"

spe_MoransI <- readRDS(here("outputs", "results", "simulations",sim_dir, "spe_MoransI.rds"))
i=1
spe_nnSVG <- readRDS(here("outputs", "simulations", sim_dir, paste0("spe_nnSVG_", i, ".rds")))
spe_SPARKX <- readRDS(here("outputs", "results", "simulations",sim_dir, "spe_SPARKX.rds"))
spe_SpaGFT <- readRDS(here("outputs", "results", "simulations",sim_dir, "spe_SpaGFT.rds"))
spe_SpatialDE2 <- readRDS(here("outputs", "results", "simulations",sim_dir, "spe_SpatialDE2.rds"))
n_genes <- dim(spe_HVG)[1]

df_HVG <- data.frame(
  rank = rowData(spe_HVG)$HVG_rank,
  mean = rowData(spe_HVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column() 

HVG_ridge_combined <- df_HVG %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = rank, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "rank",
    title = "HVG"
  ) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw()

df_MoransI <- data.frame(
  rank = rowData(spe_MoransI)$MoransI_rank,
  #moransI did not calculate mean
  mean = rowData(spe_HVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

MoransI_ridge_combined <- df_MoransI %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = rank, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "rank",
    title = "MoransI"
  ) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw()

df_nnSVG <- data.frame(
  rank = rowData(spe_nnSVG)$rank,
  mean = rowData(spe_nnSVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

nnSVG_ridge_combined <- df_nnSVG %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = rank, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "rank",
    title = "nnSVG"
  ) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw()

df_SPARKX <- data.frame(
  rank = rowData(spe_SPARKX)$SPARKX_rank,
  #sparkX did not calculate mean
  mean = rowData(spe_HVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

SPARKX_ridge_combined <- df_SPARKX %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = rank, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "rank",
    title = "SPARKX"
  ) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw()

df_SpaGFT <- data.frame(
  rank = rowData(spe_SpaGFT)$SpaGFT_rank,
  #SpaGFT did not calculate mean
  mean = rowData(spe_HVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

SpaGFT_ridge_combined <- df_SpaGFT %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = rank, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "rank",
    title = "SpaGFT"
  ) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw()

df_SpatialDE2 <- data.frame(
  rank = rowData(spe_SpatialDE2)$SpatialDE2_rank,
  #SpaGFT did not calculate mean
  mean = rowData(spe_HVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

SpatialDE2_ridge_combined <- df_SpatialDE2 %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = rank, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "rank",
    title = "SpatialDE2"
  ) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw()

ggsave(here("plots", "main", "simulations_comparing_SVG_methods_separated.png"),
       wrap_plots(
          HVG_ridge_combined,
          MoransI_ridge_combined,
          nnSVG_ridge_combined,
          SPARKX_ridge_combined,
          SpaGFT_ridge_combined,
          SpatialDE2_ridge_combined,
          ncol = 3
), width = 21, height = 20)
