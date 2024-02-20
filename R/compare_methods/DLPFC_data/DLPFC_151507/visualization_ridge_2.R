library(ggplot2)
library(SpatialExperiment)
library(patchwork)
library(GGally)
library(dplyr)
library(ggridges)

spe_HVG <- readRDS("spe_HVG.rds")
spe_MoransI <- readRDS("spe_MoransI.rds")
spe_nnSVG <- readRDS("spe_nnSVG.rds")
spe_SPARKX <- readRDS("spe_SPARKX.rds")

SpaGFT_results <- read.table("/users/kshah/MeanVarBias/R/compare_methods/DLPFC_data/DLPFC_151507/SpaGFT_results.csv",
                              header = T, sep = ",", row.names = 1)
SpatialDE2_results <- read.table("/users/kshah/MeanVarBias/R/compare_methods/DLPFC_data/DLPFC_151507/SpatialDE2_results.csv",
                                 header = T, sep = ",", row.names = 1)
n_genes <- dim(spe_nnSVG)[1]

#ridge plots separated by noise and signal for HVG
frac <- round(dim(spe_HVG)[1]*0.1*0.1)

df <- data.frame(
  rank = rowData(spe_HVG)$HVG_rank,
  mean = rowData(spe_HVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column() 

HVG_ridge_combined <- df %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = rank, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "rank",
    title = "HVG combined"
  ) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw()

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

HVG_ridge_signal <- ggplot(df, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3,
                       stat = "binline", bins = 75, draw_baseline = FALSE) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile - mean of logcounts",
    x = "rank",
    title = "HVG signal"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

HVG_ridge <- wrap_plots(HVG_ridge_combined, HVG_ridge_signal)

#ridge plots separated by noise and signal for MoransI
frac <- round(dim(spe_MoransI)[1]*0.1*0.1)

df <- data.frame(
  rank = rowData(spe_MoransI)$MoransI_rank,
  #moransI did not calculate mean
  mean = rowData(spe_HVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

MoransI_ridge_combined <- df %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = rank, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "rank",
    title = "MoransI combined"
  ) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw()

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

MoransI_ridge_signal <- ggplot(df, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3,
                       stat = "binline", bins = 75, draw_baseline = FALSE) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile - mean of logcounts",
    x = "rank",
    title = "MoransI signal"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

MoransI_ridge <- wrap_plots(MoransI_ridge_combined, MoransI_ridge_signal)

#ridge plots separated by noise and signal for nnSVG
frac <- round(dim(spe_nnSVG)[1]*0.1*0.1)

df <- data.frame(
  rank = rowData(spe_nnSVG)$rank,
  mean = rowData(spe_nnSVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()


nnSVG_ridge_combined <- df %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = rank, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "rank",
    title = "nnSVG combined"
  ) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw()

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

nnSVG_ridge_signal <- ggplot(df, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3,
                       stat = "binline", bins = 75, draw_baseline = FALSE) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile - mean of logcounts",
    x = "rank",
    title = "nnSVG signal"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

nnSVG_ridge <- wrap_plots(nnSVG_ridge_combined, nnSVG_ridge_signal)

#ridge plots separated by noise and signal for SPARKX
frac <- round(dim(spe_SPARKX)[1]*0.1*0.1)

df <- data.frame(
  rank = rowData(spe_SPARKX)$SPARKX_rank,
  #sparkX did not calculate mean
  mean = rowData(spe_HVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()


SPARKX_ridge_combined <- df %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = rank, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "rank",
    title = "SPARKX combined"
  ) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw()

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

SPARKX_ridge_signal <- ggplot(df, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3,
                       stat = "binline", bins = 75, draw_baseline = FALSE) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile - mean of logcounts",
    x = "rank",
    title = "SPARKX signal"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

SPARKX_ridge <- wrap_plots(SPARKX_ridge_combined, SPARKX_ridge_signal)

#ridge plots separated by noise and signal for SpaGFT
frac <- round(dim(spe_HVG)[1]*0.1*0.1)

df_merge <- merge(rowData(spe_nnSVG), SpaGFT_results, by=0)
df_merge$SpaGFT_rank <- rank(df_merge$qvalue, ties.method = "first")
row.names(df_merge) <- df_merge$Row.names

df <- data.frame(
  rank = df_merge$SpaGFT_rank,
  mean = df_merge$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

SpaGFT_ridge_combined <- df %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = rank, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "rank",
    title = "SpaGFT combined"
  ) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw()

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

SpaGFT_ridge_signal <- ggplot(df, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3,
                       stat = "binline", bins = 75, draw_baseline = FALSE) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile - mean of logcounts",
    x = "rank",
    title = "SpaGFT signal"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

SpaGFT_ridge <- wrap_plots(SpaGFT_ridge_combined, SpaGFT_ridge_signal)

#ridge plots separated by noise and signal for SpatialDE2
frac <- round(dim(spe_HVG)[1]*0.1*0.1)

df_merge <- merge(rowData(spe_nnSVG), SpatialDE2_results, by=0)

# rank largest fraction of spatial variances on top
df_merge$SpatialDE2_rank <- rank(-df_merge$FSV, ties.method = "first")
row.names(df_merge) <- df_merge$Row.names

df <- data.frame(
  rank = df_merge$SpatialDE2_rank,
  mean = df_merge$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

SpatialDE2_ridge_combined <- df %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = rank, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "rank",
    title = "SpatialDE2 combined"
  ) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw()

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

SpatialDE2_ridge_signal <- ggplot(df, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3,
                       stat = "binline", bins = 75, draw_baseline = FALSE) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile - mean of logcounts",
    x = "rank",
    title = "SpatialDE2 signal"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

SpatialDE2_ridge <- wrap_plots(SpatialDE2_ridge_combined, SpatialDE2_ridge_signal)


pdf("visualization_ridge_2.pdf", width = 21, height = 20)
HVG_ridge
MoransI_ridge
nnSVG_ridge
SPARKX_ridge
dev.off()
