library(ggplot2)
library(SpatialExperiment)
library(patchwork)
library(GGally)
library(dplyr)
library(ggridges)
library(here)

load("MeanVarBias/R/simulations/all_SVGs_300/simulation_weighted_SPARKX.Rdata")

unweighted_rank <- rank(sparkX_unweighted$res_mtest["adjustedPval"])
weighted_rank <- rank(sparkX_weighted$res_mtest["adjustedPval"])

spe <- readRDS(file = "MeanVarBias/R/simulations/all_SVGs_300/spe_simulation.rds")

unweighted_mean <- rowMeans(counts(spe))

#overlay unweighted and weighted ridge plots
df_unw <- data.frame(
  rank = unweighted_rank,
  mean = unweighted_mean,
  method = rep("unw", 300) 
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_w <- data.frame(
  rank = weighted_rank,
  mean = unweighted_mean,
  method = rep("w", 300) 
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df <- rbind(df_unw, df_w) %>% 
  mutate(quantile = as.factor(quantile))

ridge_overlay <- ggplot(df, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = method), rel_min_height = 0.02, alpha = 0.3, scale = 1.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile - unw mean of counts",
    x = "rank",
    title = "Ridge plots: effect of weighting on rank"
  ) +
  scale_fill_manual(labels = c("weighted", "unweighted"), values = c("red", "blue")) +
  coord_cartesian(xlim = c(1, 300)) +
  theme_bw()

ggsave(here("MeanVarBias", "R", "simulations", "all_SVGs_300", "ridge_overlay_SPARKX.png"), ridge_overlay)

#ridge plots separated by noise and signal for unweighted and weighted
frac <- round(dim(spe)[1]*0.1)*0.1

df_unw_signal <- df_unw %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "signal")

indices <- as.integer(df_unw_signal$rowname)

df_unw_background <- df_unw[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "background")

df <- rbind(df_unw_signal, df_unw_background)

rank_separated_unw <- ggplot(df, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile - unw mean of logcounts",
    x = "rank",
    title = "Signal unweighted"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, 300)) +
  theme_bw() 

df_w_signal <- df_w %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "signal")

indices <- as.integer(df_w_signal$rowname)

df_w_background <- df_w[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "background")

df <- rbind(df_w_signal, df_w_background)

rank_separated_w <- ggplot(df, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "decile - unw mean of logcounts",
    x = "rank",
    title = "Signal weighted"
  )  +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, 300)) +
  theme_bw() 

ridge_signal <- wrap_plots(rank_separated_unw, rank_separated_w, nrow=1, guides = "collect") 

ggsave(here("MeanVarBias", "R", "simulations", "all_SVGs_300", "ridge_signal_SPARKX.png"), ridge_signal)
