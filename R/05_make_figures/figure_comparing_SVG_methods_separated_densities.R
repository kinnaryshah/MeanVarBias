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
spe_SMASH <- readRDS(here("outputs", "results", "spe_humanDLPFC_SMASH.rds"))
spe_HEARTSVG <- readRDS(here("outputs", "results", "spe_humanDLPFC_HEARTSVG.rds"))

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

df_MoransI$quantile <- factor(df_MoransI$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

MoransI_ridge <- ggplot(df_MoransI, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "DLPFC Deciles",
    x = "Rank",
    title = "Moran's I"
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
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

df_nnSVG$quantile <- factor(df_nnSVG$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

nnSVG_ridge <- ggplot(df_nnSVG, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "DLPFC Deciles",
    x = "Rank",
    title = "nnSVG"
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
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

df_SPARKX$quantile <- factor(df_SPARKX$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

SPARKX_ridge <- ggplot(df_SPARKX, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "DLPFC Deciles",
    x = "Rank",
    title = "SPARK-X"
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
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

df_SpaGFT$quantile <- factor(df_SpaGFT$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

SpaGFT_ridge <- ggplot(df_SpaGFT, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "DLPFC Deciles",
    x = "Rank",
    title = "SpaGFT"
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
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

df_SpatialDE2$quantile <- factor(df_SpatialDE2$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

SpatialDE2_ridge <- ggplot(df_SpatialDE2, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "DLPFC Deciles",
    x = "Rank",
    title = "SpatialDE2"
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  My_Theme

df_SMASH <- data.frame(
  rank = rowData(spe_SMASH)$SMASH_rank,
  #SpaGFT did not calculate mean
  mean = rowData(spe_nnSVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_SMASH_signal <- df_SMASH %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "Signal")

indices <- as.integer(df_SMASH_signal$rowname)

df_SMASH_background <- df_SMASH[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_SMASH <- rbind(df_SMASH_signal, df_SMASH_background)

df_SMASH$quantile <- factor(df_SMASH$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

SMASH_ridge <- ggplot(df_SMASH, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "DLPFC Deciles",
    x = "Rank",
    title = "SMASH"
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  My_Theme

df_HEARTSVG <- data.frame(
  rank = rowData(spe_HEARTSVG)$HEARTSVG_rank,
  #SpaGFT did not calculate mean
  mean = rowData(spe_nnSVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_HEARTSVG_signal <- df_HEARTSVG %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "Signal")

indices <- as.integer(df_HEARTSVG_signal$rowname)

df_HEARTSVG_background <- df_HEARTSVG[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_HEARTSVG <- rbind(df_HEARTSVG_signal, df_HEARTSVG_background)

df_HEARTSVG$quantile <- factor(df_HEARTSVG$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

HEARTSVG_ridge <- ggplot(df_HEARTSVG, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "DLPFC Deciles",
    x = "Rank",
    title = "HEARTSVG"
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  My_Theme


DLPFC_plotlist <- list(
  MoransI_ridge,
  nnSVG_ridge,
  SPARKX_ridge,
  SpaGFT_ridge,
  SpatialDE2_ridge,
  SMASH_ridge,
  HEARTSVG_ridge
)

spe_MoransI <- readRDS(here("outputs", "results", "spe_humanOvarian_MoransI.rds"))
spe_nnSVG <- readRDS(here("outputs", "results", "spe_humanOvarian_nnSVG.rds"))
spe_SPARKX <- readRDS(here("outputs", "results", "spe_humanOvarian_SPARKX.rds"))
spe_SpaGFT <- readRDS(here("outputs", "results", "spe_humanOvarian_SpaGFT.rds"))
spe_SpatialDE2 <- readRDS(here("outputs", "results", "spe_humanOvarian_SpatialDE2.rds"))

# something wrong with files below
# need to run run_SMASH_Ovarian.R and save spe 
spe_SMASH <- readRDS(here("outputs", "results", "spe_humanOvarian_SMASH.rds"))
spe_HEARTSVG <- readRDS(here("outputs", "results", "spe_humanOvarian_HEARTSVG.rds"))

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

df_MoransI$quantile <- factor(df_MoransI$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

Ovarian_MoransI_ridge <- ggplot(df_MoransI, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Ovarian Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
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

df_nnSVG$quantile <- factor(df_nnSVG$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

Ovarian_nnSVG_ridge <- ggplot(df_nnSVG, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Ovarian Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
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

df_SPARKX$quantile <- factor(df_SPARKX$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

Ovarian_SPARKX_ridge <- ggplot(df_SPARKX, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Ovarian Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
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

df_SpaGFT$quantile <- factor(df_SpaGFT$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

Ovarian_SpaGFT_ridge <- ggplot(df_SpaGFT, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Ovarian Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
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

df_SpatialDE2$quantile <- factor(df_SpatialDE2$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

Ovarian_SpatialDE2_ridge <- ggplot(df_SpatialDE2, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Ovarian Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  My_Theme

df_SMASH <- data.frame(
  rank = rowData(spe_SMASH)$SMASH_rank,
  #SpaGFT did not calculate mean
  mean = rowData(spe_nnSVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_SMASH_signal <- df_SMASH %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "Signal")

indices <- as.integer(df_SMASH_signal$rowname)

df_SMASH_background <- df_SMASH[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_SMASH <- rbind(df_SMASH_signal, df_SMASH_background)

df_SMASH$quantile <- factor(df_SMASH$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

Ovarian_SMASH_ridge <- ggplot(df_SMASH, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Ovarian Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  My_Theme

df_HEARTSVG <- data.frame(
  rank = rowData(spe_HEARTSVG)$HEARTSVG_rank,
  #SpaGFT did not calculate mean
  mean = rowData(spe_nnSVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_HEARTSVG_signal <- df_HEARTSVG %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "Signal")

indices <- as.integer(df_HEARTSVG_signal$rowname)

df_HEARTSVG_background <- df_HEARTSVG[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_HEARTSVG <- rbind(df_HEARTSVG_signal, df_HEARTSVG_background)

df_HEARTSVG$quantile <- factor(df_HEARTSVG$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

Ovarian_HEARTSVG_ridge <- ggplot(df_HEARTSVG, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Ovarian Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  My_Theme


Ovarian_plotlist <- list(
  Ovarian_MoransI_ridge,
  Ovarian_nnSVG_ridge,
  Ovarian_SPARKX_ridge,
  Ovarian_SpaGFT_ridge,
  Ovarian_SpatialDE2_ridge,
  Ovarian_SMASH_ridge,
  Ovarian_HEARTSVG_ridge
)

spe_MoransI <- readRDS(here("outputs", "results", "spe_humanLobularBreast_MoransI.rds"))
spe_nnSVG <- readRDS(here("outputs", "results", "spe_humanLobularBreast_nnSVG.rds"))
spe_SPARKX <- readRDS(here("outputs", "results", "spe_humanLobularBreast_SPARKX.rds"))
spe_SpaGFT <- readRDS(here("outputs", "results", "spe_humanLobularBreast_SpaGFT.rds"))
spe_SpatialDE2 <- readRDS(here("outputs", "results", "spe_humanLobularBreast_SpatialDE2.rds"))
spe_SMASH <- readRDS(here("outputs", "results", "spe_humanLobularBreast_SMASH.rds"))
spe_HEARTSVG <- readRDS(here("outputs", "results", "spe_humanLobularBreast_HEARTSVG.rds"))

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

df_MoransI$quantile <- factor(df_MoransI$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

LobularBreast_MoransI_ridge <- ggplot(df_MoransI, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Lobular Breast Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
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

df_nnSVG$quantile <- factor(df_nnSVG$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

LobularBreast_nnSVG_ridge <- ggplot(df_nnSVG, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Lobular Breast Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
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

df_SPARKX$quantile <- factor(df_SPARKX$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

LobularBreast_SPARKX_ridge <- ggplot(df_SPARKX, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Lobular Breast Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
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

df_SpaGFT$quantile <- factor(df_SpaGFT$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

LobularBreast_SpaGFT_ridge <- ggplot(df_SpaGFT, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Lobular Breast Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
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

df_SpatialDE2$quantile <- factor(df_SpatialDE2$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

LobularBreast_SpatialDE2_ridge <- ggplot(df_SpatialDE2, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Lobular Breast Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  My_Theme

df_SMASH <- data.frame(
  rank = rowData(spe_SMASH)$SMASH_rank,
  #SpaGFT did not calculate mean
  mean = rowData(spe_nnSVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_SMASH_signal <- df_SMASH %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "Signal")

indices <- as.integer(df_SMASH_signal$rowname)

df_SMASH_background <- df_SMASH[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_SMASH <- rbind(df_SMASH_signal, df_SMASH_background)

df_SMASH$quantile <- factor(df_SMASH$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

LobularBreast_SMASH_ridge <- ggplot(df_SMASH, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Lobular Breast Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  My_Theme

df_HEARTSVG <- data.frame(
  rank = rowData(spe_HEARTSVG)$HEARTSVG_rank,
  #SpaGFT did not calculate mean
  mean = rowData(spe_nnSVG)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_HEARTSVG_signal <- df_HEARTSVG %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "Signal")

indices <- as.integer(df_HEARTSVG_signal$rowname)

df_HEARTSVG_background <- df_HEARTSVG[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_HEARTSVG <- rbind(df_HEARTSVG_signal, df_HEARTSVG_background)

df_HEARTSVG$quantile <- factor(df_HEARTSVG$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

LobularBreast_HEARTSVG_ridge <- ggplot(df_HEARTSVG, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Lobular Breast Deciles",
    x = "Rank",
    title = ""
  ) +
  guides(fill=guide_legend(title="Group")) +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text = element_text(hjust = 0),
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_rect(colour="black", fill="white")
  ) +
  My_Theme

LobularBreast_plotlist <- list(
  LobularBreast_MoransI_ridge,
  LobularBreast_nnSVG_ridge,
  LobularBreast_SPARKX_ridge,
  LobularBreast_SpaGFT_ridge,
  LobularBreast_SpatialDE2_ridge,
  LobularBreast_SMASH_ridge,
  LobularBreast_HEARTSVG_ridge
)

ggsave(here("plots", "main", "comparing_SVG_methods_separated.png"),
       wrap_plots(
         DLPFC_plotlist[[1]],
         DLPFC_plotlist[[2]],
         DLPFC_plotlist[[3]],
         DLPFC_plotlist[[4]],
         DLPFC_plotlist[[5]],
         DLPFC_plotlist[[6]],
         DLPFC_plotlist[[7]],
         Ovarian_plotlist[[1]],
         Ovarian_plotlist[[2]],
         Ovarian_plotlist[[3]],
         Ovarian_plotlist[[4]],
         Ovarian_plotlist[[5]],
         Ovarian_plotlist[[6]],
         Ovarian_plotlist[[7]],
         LobularBreast_plotlist[[1]],
         LobularBreast_plotlist[[2]],
         LobularBreast_plotlist[[3]],
         LobularBreast_plotlist[[4]],
         LobularBreast_plotlist[[5]],
         LobularBreast_plotlist[[6]],
         LobularBreast_plotlist[[7]],
         nrow = 3,
         guides = "collect",
         axis_titles = "collect"
       ) + plot_annotation(tag_levels = 'A') & 
         theme(plot.tag = element_text(size = 25)),
       width = 25, height = 25)
