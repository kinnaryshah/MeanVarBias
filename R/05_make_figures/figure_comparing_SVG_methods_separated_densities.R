library(ggplot2)
library(SpatialExperiment)
library(patchwork)
library(GGally)
library(dplyr)
library(ggridges)
library(here)


My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  legend.title = element_text(size=20),
  legend.text = element_text(size=20),
  title = element_text(size = 20))


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
  mutate(grp = "Signal")

indices <- as.integer(df_MoransI_signal$rowname)

df_MoransI_background <- df_MoransI[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_MoransI <- rbind(df_MoransI_signal, df_MoransI_background)

MoransI_ridge <- ggplot(df_MoransI, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "DLPFC Deciles",
    x = "Rank",
    title = "Moran's I"
  ) +
  guides(fill=guide_legend(title="Group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() +
  My_Theme

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
  mutate(grp = "Signal")

indices <- as.integer(df_nnSVG_signal$rowname)

df_nnSVG_background <- df_nnSVG[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_nnSVG <- rbind(df_nnSVG_signal, df_nnSVG_background)

nnSVG_ridge <- ggplot(df_nnSVG, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "DLPFC Deciles",
    x = "Rank",
    title = "nnSVG"
  ) +
  guides(fill=guide_legend(title="Group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() +
  My_Theme

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
  mutate(grp = "Signal")

indices <- as.integer(df_SPARKX_signal$rowname)

df_SPARKX_background <- df_SPARKX[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_SPARKX <- rbind(df_SPARKX_signal, df_SPARKX_background)

SPARKX_ridge <- ggplot(df_SPARKX, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "DLPFC Deciles",
    x = "Rank",
    title = "SPARK-X"
  ) +
  guides(fill=guide_legend(title="Group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() +
  My_Theme

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
  mutate(grp = "Signal")

indices <- as.integer(df_SpaGFT_signal$rowname)

df_SpaGFT_background <- df_SpaGFT[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_SpaGFT <- rbind(df_SpaGFT_signal, df_SpaGFT_background)

SpaGFT_ridge <- ggplot(df_SpaGFT, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "DLPFC Deciles",
    x = "Rank",
    title = "SpaGFT"
  ) +
  guides(fill=guide_legend(title="Group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() +
  My_Theme

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
  mutate(grp = "Signal")

indices <- as.integer(df_SpatialDE2_signal$rowname)

df_SpatialDE2_background <- df_SpatialDE2[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_SpatialDE2 <- rbind(df_SpatialDE2_signal, df_SpatialDE2_background)

SpatialDE2_ridge <- ggplot(df_SpatialDE2, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "DLPFC Deciles",
    x = "Rank",
    title = "SpatialDE2"
  ) +
  guides(fill=guide_legend(title="Group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() +
  My_Theme

DLPFC_plotlist <- list(
  MoransI_ridge,
  nnSVG_ridge,
  SPARKX_ridge,
  SpaGFT_ridge,
  SpatialDE2_ridge
)

spe_MoransI <- readRDS(here("outputs", "results", "spe_humanOvarian_MoransI.rds"))
spe_nnSVG <- readRDS(here("outputs", "results", "spe_humanOvarian_nnSVG.rds"))
spe_SPARKX <- readRDS(here("outputs", "results", "spe_humanOvarian_SPARKX.rds"))
spe_SpaGFT <- readRDS(here("outputs", "results", "spe_humanOvarian_SpaGFT.rds"))
spe_SpatialDE2 <- readRDS(here("outputs", "results", "spe_humanOvarian_SpatialDE2.rds"))
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
  mutate(grp = "Signal")

indices <- as.integer(df_MoransI_signal$rowname)

df_MoransI_background <- df_MoransI[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_MoransI <- rbind(df_MoransI_signal, df_MoransI_background)

Ovarian_MoransI_ridge <- ggplot(df_MoransI, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "Ovarian Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() +
  My_Theme

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
  mutate(grp = "Signal")

indices <- as.integer(df_nnSVG_signal$rowname)

df_nnSVG_background <- df_nnSVG[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_nnSVG <- rbind(df_nnSVG_signal, df_nnSVG_background)

Ovarian_nnSVG_ridge <- ggplot(df_nnSVG, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "Ovarian Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() +
  My_Theme

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
  mutate(grp = "Signal")

indices <- as.integer(df_SPARKX_signal$rowname)

df_SPARKX_background <- df_SPARKX[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_SPARKX <- rbind(df_SPARKX_signal, df_SPARKX_background)

Ovarian_SPARKX_ridge <- ggplot(df_SPARKX, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "Ovarian Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() +
  My_Theme

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
  mutate(grp = "Signal")

indices <- as.integer(df_SpaGFT_signal$rowname)

df_SpaGFT_background <- df_SpaGFT[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_SpaGFT <- rbind(df_SpaGFT_signal, df_SpaGFT_background)

Ovarian_SpaGFT_ridge <- ggplot(df_SpaGFT, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "Ovarian Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() +
  My_Theme

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
  mutate(grp = "Signal")

indices <- as.integer(df_SpatialDE2_signal$rowname)

df_SpatialDE2_background <- df_SpatialDE2[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_SpatialDE2 <- rbind(df_SpatialDE2_signal, df_SpatialDE2_background)

Ovarian_SpatialDE2_ridge <- ggplot(df_SpatialDE2, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "Ovarian Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() +
  My_Theme

Ovarian_plotlist <- list(
  Ovarian_MoransI_ridge,
  Ovarian_nnSVG_ridge,
  Ovarian_SPARKX_ridge,
  Ovarian_SpaGFT_ridge,
  Ovarian_SpatialDE2_ridge
)

spe_MoransI <- readRDS(here("outputs", "results", "spe_humanLobularBreast_MoransI.rds"))
spe_nnSVG <- readRDS(here("outputs", "results", "spe_humanLobularBreast_nnSVG.rds"))
spe_SPARKX <- readRDS(here("outputs", "results", "spe_humanLobularBreast_SPARKX.rds"))
spe_SpaGFT <- readRDS(here("outputs", "results", "spe_humanLobularBreast_SpaGFT.rds"))
spe_SpatialDE2 <- readRDS(here("outputs", "results", "spe_humanLobularBreast_SpatialDE2.rds"))
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
  mutate(grp = "Signal")

indices <- as.integer(df_MoransI_signal$rowname)

df_MoransI_background <- df_MoransI[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_MoransI <- rbind(df_MoransI_signal, df_MoransI_background)

LobularBreast_MoransI_ridge <- ggplot(df_MoransI, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "Lobular Breast Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() +
  My_Theme

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
  mutate(grp = "Signal")

indices <- as.integer(df_nnSVG_signal$rowname)

df_nnSVG_background <- df_nnSVG[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_nnSVG <- rbind(df_nnSVG_signal, df_nnSVG_background)

LobularBreast_nnSVG_ridge <- ggplot(df_nnSVG, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "Lobular Breast Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() +
  My_Theme

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
  mutate(grp = "Signal")

indices <- as.integer(df_SPARKX_signal$rowname)

df_SPARKX_background <- df_SPARKX[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_SPARKX <- rbind(df_SPARKX_signal, df_SPARKX_background)

LobularBreast_SPARKX_ridge <- ggplot(df_SPARKX, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "Lobular Breast Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() +
  My_Theme

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
  mutate(grp = "Signal")

indices <- as.integer(df_SpaGFT_signal$rowname)

df_SpaGFT_background <- df_SpaGFT[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_SpaGFT <- rbind(df_SpaGFT_signal, df_SpaGFT_background)

LobularBreast_SpaGFT_ridge <- ggplot(df_SpaGFT, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "Lobular Breast Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() +
  My_Theme

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
  mutate(grp = "Signal")

indices <- as.integer(df_SpatialDE2_signal$rowname)

df_SpatialDE2_background <- df_SpatialDE2[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_SpatialDE2 <- rbind(df_SpatialDE2_signal, df_SpatialDE2_background)

LobularBreast_SpatialDE2_ridge <- ggplot(df_SpatialDE2, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "Lobular Breast Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() +
  My_Theme

LobularBreast_plotlist <- list(
  LobularBreast_MoransI_ridge,
  LobularBreast_nnSVG_ridge,
  LobularBreast_SPARKX_ridge,
  LobularBreast_SpaGFT_ridge,
  LobularBreast_SpatialDE2_ridge
)


ggsave(here("plots", "main", "comparing_SVG_methods_separated.png"),
       wrap_plots(
         DLPFC_plotlist[[1]],
         DLPFC_plotlist[[2]],
         DLPFC_plotlist[[3]],
         DLPFC_plotlist[[4]],
         DLPFC_plotlist[[5]],
         Ovarian_plotlist[[1]],
         Ovarian_plotlist[[2]],
         Ovarian_plotlist[[3]],
         Ovarian_plotlist[[4]],
         Ovarian_plotlist[[5]],
         LobularBreast_plotlist[[1]],
         LobularBreast_plotlist[[2]],
         LobularBreast_plotlist[[3]],
         LobularBreast_plotlist[[4]],
         LobularBreast_plotlist[[5]],
         nrow = 3,
         guides = "collect",
         axis_titles = "collect"
       ) + plot_annotation(tag_levels = 'A') & 
         theme(plot.tag = element_text(size = 25)),
       width = 21, height = 30)
