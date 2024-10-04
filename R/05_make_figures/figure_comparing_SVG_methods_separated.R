library(ggplot2)
library(SpatialExperiment)
library(patchwork)
library(GGally)
library(dplyr)
library(ggridges)
library(here)

spe_MoransI <- readRDS(here("outputs", "results", "spe_humanDLPFC_MoransI.rds"))
spe_nnSVG <- readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG.rds"))
spe_SPARKX <- readRDS(here("outputs", "results", "spe_humanDLPFC_SPARKX.rds"))
spe_SpaGFT <- readRDS(here("outputs", "results", "spe_humanDLPFC_SpaGFT.rds"))
spe_SpatialDE2 <- readRDS(here("outputs", "results", "spe_humanDLPFC_SpatialDE2.rds"))
n_genes <- dim(spe_MoransI)[1]
frac <- round(n_genes*0.1*0.1)

df_MoransI <- data.frame(
  rank = rowData(spe_MoransI)$MoransI_rank,
  #moransI did not calculate mean
  mean = rowData(spe_nnSVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_MoransI_signal <- df_MoransI %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "signal")

indices <- as.integer(df_MoransI_signal$rowname)

df_MoransI_background <- df_MoransI[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "background")

df_MoransI <- rbind(df_MoransI_signal, df_MoransI_background)

MoransI_ridge <- ggplot(df_MoransI, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile = mean of logcounts",
    x = "rank",
    title = "Moran's I"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

df_nnSVG <- data.frame(
  rank = rowData(spe_nnSVG)$rank,
  mean = rowData(spe_nnSVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_nnSVG_signal <- df_nnSVG %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "signal")

indices <- as.integer(df_nnSVG_signal$rowname)

df_nnSVG_background <- df_nnSVG[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "background")

df_nnSVG <- rbind(df_nnSVG_signal, df_nnSVG_background)

nnSVG_ridge <- ggplot(df_nnSVG, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile = mean of logcounts",
    x = "rank",
    title = "nnSVG"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

df_SPARKX <- data.frame(
  rank = rowData(spe_SPARKX)$SPARKX_rank,
  #sparkX did not calculate mean
  mean = rowData(spe_nnSVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_SPARKX_signal <- df_SPARKX %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "signal")

indices <- as.integer(df_SPARKX_signal$rowname)

df_SPARKX_background <- df_SPARKX[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "background")

df_SPARKX <- rbind(df_SPARKX_signal, df_SPARKX_background)

SPARKX_ridge <- ggplot(df_SPARKX, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile = mean of logcounts",
    x = "rank",
    title = "SPARK-X"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

df_SpaGFT <- data.frame(
  rank = rowData(spe_SpaGFT)$SpaGFT_rank,
  #SpaGFT did not calculate mean
  mean = rowData(spe_nnSVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_SpaGFT_signal <- df_SpaGFT %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "signal")

indices <- as.integer(df_SpaGFT_signal$rowname)

df_SpaGFT_background <- df_SpaGFT[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "background")

df_SpaGFT <- rbind(df_SpaGFT_signal, df_SpaGFT_background)

SpaGFT_ridge <- ggplot(df_SpaGFT, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile = mean of logcounts",
    x = "rank",
    title = "SpaGFT"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

df_SpatialDE2 <- data.frame(
  rank = rowData(spe_SpatialDE2)$SpatialDE2_rank,
  #SpaGFT did not calculate mean
  mean = rowData(spe_nnSVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_SpatialDE2_signal <- df_SpatialDE2 %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "signal")

indices <- as.integer(df_SpatialDE2_signal$rowname)

df_SpatialDE2_background <- df_SpatialDE2[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "background")

df_SpatialDE2 <- rbind(df_SpatialDE2_signal, df_SpatialDE2_background)

SpatialDE2_ridge <- ggplot(df_SpatialDE2, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile = mean of logcounts",
    x = "rank",
    title = "SpatialDE2"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

ggsave(here("plots", "main", "simulations_comparing_SVG_methods_separated.png"),
       wrap_plots(
         MoransI_ridge,
         nnSVG_ridge,
         SPARKX_ridge,
         SpaGFT_ridge,
         SpatialDE2_ridge,
         nrow = 1,
         guides = "collect",
         axis_titles = "collect"
       ), width = 21, height = 10)
