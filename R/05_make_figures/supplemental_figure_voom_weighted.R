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

spe_DLPFC <- readRDS(here("outputs", "results", "spe_humanDLPFC_voom_weighted_nnSVG.rds"))
spe_DLPFC_unweighted <- readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG.rds"))

n_genes <- dim(spe_DLPFC)[1]
frac <- round(n_genes*0.1*0.1)

df_DLPFC <- data.frame(
  rank = rowData(spe_DLPFC)$weighted_rank,
  mean = rowData(spe_DLPFC_unweighted)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_DLPFC_signal <- df_DLPFC %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "Signal")

indices <- as.integer(df_DLPFC_signal$rowname)

df_DLPFC_background <- df_DLPFC[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_DLPFC <- rbind(df_DLPFC_signal, df_DLPFC_background)

df_DLPFC$quantile <- factor(df_DLPFC$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

DLPFC_ridge <- ggplot(df_DLPFC, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "DLPFC Deciles",
    x = "Rank",
    title = "voom nnSVG"
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


spe_DLPFC <- readRDS(here("outputs", "results", "spe_humanDLPFC_weighted_nnSVG.rds"))
n_genes <- dim(spe_DLPFC)[1]
frac <- round(n_genes*0.1*0.1)

df_DLPFC <- data.frame(
  rank = rowData(spe_DLPFC)$weighted_rank,
  mean = rowData(spe_DLPFC_unweighted)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_DLPFC_signal <- df_DLPFC %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "Signal")

indices <- as.integer(df_DLPFC_signal$rowname)

df_DLPFC_background <- df_DLPFC[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_DLPFC <- rbind(df_DLPFC_signal, df_DLPFC_background)

df_DLPFC$quantile <- factor(df_DLPFC$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

DLPFC_ridge_spoon <- ggplot(df_DLPFC, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "DLPFC Deciles",
    x = "Rank",
    title = "spoon nnSVG"
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


file_dir <- here("outputs", "simulations", "reps_968_100_50per_1000_0.2_to_3_0.5_to_9")
spe_sim <- readRDS(paste0(file_dir, "/spe_voom_weighted_nnSVG_1.rds"))
spe_sim_unweighted <- readRDS(paste0(file_dir, "/spe_nnSVG_1.rds"))

n_genes <- dim(spe_sim)[1]
frac <- round(n_genes*0.1*0.1)

df_sim <- data.frame(
  rank = rowData(spe_sim)$weighted_rank,
  mean = rowData(spe_sim_unweighted)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_sim_signal <- df_sim %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "Signal")

indices <- as.integer(df_sim_signal$rowname)

df_sim_background <- df_sim[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_sim <- rbind(df_sim_signal, df_sim_background)

df_sim$quantile <- factor(df_sim$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

sim_ridge <- ggplot(df_sim, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Simulation Deciles",
    x = "Rank",
    title = "voom nnSVG"
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

file_dir <- here("outputs", "simulations", "reps_968_100_50per_1000_0.2_to_3_0.5_to_9")
spe_sim <- readRDS(paste0(file_dir, "/spe_weighted_nnSVG_1.rds"))

n_genes <- dim(spe_sim)[1]
frac <- round(n_genes*0.1*0.1)

df_sim <- data.frame(
  rank = rowData(spe_sim)$weighted_rank,
  mean = rowData(spe_sim_unweighted)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_sim_signal <- df_sim %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "Signal")

indices <- as.integer(df_sim_signal$rowname)

df_sim_background <- df_sim[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "Background")

df_sim <- rbind(df_sim_signal, df_sim_background)

df_sim$quantile <- factor(df_sim$quantile, levels=c(10,9,8,7,6,5,4,3,2,1))

sim_ridge_spoon <- ggplot(df_sim, aes(x = rank, fill = grp)) +
  geom_density(alpha = 0.4, trim=T) +
  facet_wrap(~ quantile, scales = "free_y", ncol=1, switch = "y") +
  theme_bw() +
  labs(
    y = "Simulation Deciles",
    x = "Rank",
    title = "spoon nnSVG"
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

ggsave(here("plots", "supplementary", "comparing_voom.png"),
       wrap_plots(DLPFC_ridge, sim_ridge,
                  DLPFC_ridge_spoon, sim_ridge_spoon,
                  nrow = 2,
                  guides = "collect") + plot_annotation(tag_levels = 'A') & 
         theme(plot.tag = element_text(size = 25)),
       width = 15, height = 15)
