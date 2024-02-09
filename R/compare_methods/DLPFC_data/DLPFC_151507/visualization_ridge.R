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

n_genes <- dim(spe_nnSVG)[1]

#ridge plots separated by noise and signal for HVG
frac <- round(dim(spe_HVG)[1]*0.1*0.1)

df <- data.frame(
  rank = rowData(spe_HVG)$HVG_rank,
  mean = rowData(spe_HVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

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

HVG_ridge <- ggplot(df, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile - mean of logcounts",
    x = "rank",
    title = "HVG signal"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

#ridge plots separated by noise and signal for MoransI
frac <- round(dim(spe_MoransI)[1]*0.1*0.1)

df <- data.frame(
  rank = rowData(spe_MoransI)$MoransI_rank,
  #moransI did not calculate mean
  mean = rowData(spe_HVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

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

MoransI_ridge <- ggplot(df, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile - mean of logcounts",
    x = "rank",
    title = "MoransI signal"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

#ridge plots separated by noise and signal for nnSVG
frac <- round(dim(spe_nnSVG)[1]*0.1*0.1)

df <- data.frame(
  rank = rowData(spe_nnSVG)$rank,
  mean = rowData(spe_nnSVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

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

nnSVG_ridge <- ggplot(df, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile - mean of logcounts",
    x = "rank",
    title = "nnSVG signal"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

#ridge plots separated by noise and signal for SPARKX
frac <- round(dim(spe_SPARKX)[1]*0.1*0.1)

df <- data.frame(
  rank = rowData(spe_SPARKX)$SPARKX_rank,
  #sparkX did not calculate mean
  mean = rowData(spe_HVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

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

SPARKX_ridge <- ggplot(df, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile - mean of logcounts",
    x = "rank",
    title = "SPARKX signal"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

pdf("visualization_ridge.pdf", width = 21, height = 20)
HVG_ridge
MoransI_ridge
nnSVG_ridge
SPARKX_ridge
dev.off()
