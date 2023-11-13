# only constrained weighted

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
library(ggridges)
library(ggformula)
library(GGally)
library(ggblend)

#this is constrained weighted and unweighted
load(file = "mean_var_project/spe_V10B01_086_B1_HPC_weighted_nnSVG.rds")

pdf(file="mean_var_project/V10B01_086_B1_HPC_vis.pdf")

weighted_logcounts <- assay(spe_weighted, "weighted_logcounts")
weighted_var <- rowVars(weighted_logcounts)

df <- data.frame(LR_stat = rowData(spe_unweighted)$LR_stat,
                 mean = rowData(spe_unweighted)$mean,
                 var = rowData(spe_unweighted)$var,
                 sigma.sq = rowData(spe_unweighted)$sigma.sq,
                 tau.sq = rowData(spe_unweighted)$tau.sq,
                 prop_sv = rowData(spe_unweighted)$prop_sv)

region <- "HPC: unweighted"

# variance vs. mean
var1 <- ggplot(df, 
               aes(x = mean, y = var, color = LR_stat)) + 
  #scale_x_continuous(trans = "log10") + 
  #scale_y_continuous(trans = "log10") + 
  geom_point(size = 0.5) + 
  ylim(c(0,3)) +
  geom_smooth(method="loess", color="black", size=0.5) +
  scale_color_viridis(trans = "log10") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "mean logcounts", 
       y = "variance", 
       color = "LR stat",
       title = region) + 
  theme_bw() + 
  theme(strip.text.y.left = element_text(angle = 0))

# spatial variance vs. mean
spat_var1 <- ggplot(df, 
                    aes(x = mean, y = sigma.sq, color = LR_stat)) + 
  #scale_x_continuous(trans = "log10") + 
  #scale_y_continuous(trans = "log10") + 
  geom_point(size = 0.5) + 
  ylim(c(0,3)) +
  scale_color_viridis(trans = "log10") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "mean logcounts", 
       y = expression(paste("spatial variance ", sigma^{2})), 
       color = "LR stat",
       title = region) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.y = element_blank())

nonspat_var1 <- ggplot(df, 
                       aes(x = mean, y = tau.sq, color = LR_stat)) + 
  #scale_x_continuous(trans = "log10") + 
  #scale_y_continuous(trans = "log10") + 
  geom_point(size = 0.5) + 
  ylim(c(0,3)) +
  scale_color_viridis(trans = "log10") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "mean logcounts", 
       y = expression(paste("nonspatial variance ", tau^{2})), 
       color = "LR stat",
       title = region) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.y = element_blank())

prop1 <- ggplot(df, 
                aes(x = mean, y = prop_sv, color = LR_stat)) + 
  #scale_x_continuous(trans = "log10") + 
  geom_point(size = 0.5) + 
  scale_color_viridis(trans = "log10") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "mean logcounts", 
       y = expression(paste("prop. spatial variance ", sigma^{2}/(sigma^{2}+tau^{2}))), 
       color = "LR stat",
       title = region) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.y = element_blank())

df <- data.frame(LR_stat = rowData(spe_weighted)$weighted_LR_stat,
                 mean = rowData(spe_weighted)$weighted_mean,
                 var = weighted_var,
                 sigma.sq = rowData(spe_weighted)$weighted_sigma.sq,
                 tau.sq = rowData(spe_weighted)$weighted_tau.sq,
                 prop_sv = rowData(spe_weighted)$weighted_prop_sv)

region <- "HPC: constrained weighted"

# variance vs. mean
var <- ggplot(df, 
              aes(x = mean, y = var, color = LR_stat)) + 
  #scale_x_continuous(trans = "log10") + 
  #scale_y_continuous(trans = "log10") + 
  geom_point(size = 0.5) + 
  ylim(c(0,15)) +
  geom_smooth(method="loess", color="black", size=0.5) +
  scale_color_viridis(trans = "log10") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "mean logcounts", 
       y = "variance", 
       color = "LR stat",
       title = region) + 
  theme_bw() + 
  theme(strip.text.y.left = element_text(angle = 0))

# spatial variance vs. mean
spat_var <- ggplot(df, 
                   aes(x = mean, y = sigma.sq, color = LR_stat)) + 
  #scale_x_continuous(trans = "log10") + 
  #scale_y_continuous(trans = "log10") + 
  geom_point(size = 0.5) + 
  ylim(c(0,15)) +
  scale_color_viridis(trans = "log10") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "mean logcounts", 
       y = expression(paste("spatial variance ", sigma^{2})), 
       color = "LR stat",
       title = region) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.y = element_blank())

nonspat_var <- ggplot(df, 
                      aes(x = mean, y = tau.sq, color = LR_stat)) + 
  #scale_x_continuous(trans = "log10") + 
  #scale_y_continuous(trans = "log10") + 
  geom_point(size = 0.5) + 
  ylim(c(0,15)) +
  scale_color_viridis(trans = "log10") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "mean logcounts", 
       y = expression(paste("nonspatial variance ", tau^{2})), 
       color = "LR stat",
       title = region) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.y = element_blank())

prop <- ggplot(df, 
               aes(x = mean, y = prop_sv, color = LR_stat)) + 
  #scale_x_continuous(trans = "log10") + 
  geom_point(size = 0.5) + 
  scale_color_viridis(trans = "log10") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "mean logcounts", 
       y = expression(paste("prop. spatial variance ", sigma^{2}/(sigma^{2}+tau^{2}))), 
       color = "LR stat",
       title = region) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text.y = element_blank())

wrap_plots(var1, spat_var1, nonspat_var1, prop1,
           var, spat_var, nonspat_var, prop, guides="collect",
           ncol=4, nrow=2) +
  plot_annotation(title = "mean vs. components of variance",
                  caption = "no data cut off")

### EXPLORING CHANGING RANKS
#prop sv, sigma.sq, LR pre weighting, LR post weighting, etc.

weighted_rank <- rank(-1*rowData(spe_weighted)$weighted_LR_stat)

#new value for change in weighted rank
rank_shift <- rowData(spe_unweighted)$rank - weighted_rank

#only use positive direction of rank changes
top_20_weighted <- order(rank_shift, decreasing=TRUE)[1:20]
top_20_weighted_tf <- rep(F, length(rank_shift))
top_20_weighted_tf[top_20_weighted] <- T

rowData(spe_weighted)$gene_name[top_20_weighted_tf]

#plot shift in rank against unw mean expression 
rank1 <- data.frame(x = rowData(spe_unweighted)$mean,
                    y = rank_shift,
                    top_change = top_20_weighted_tf
) |>
  ggplot() +
  geom_point(aes(x = x, y = y, color = top_change), size=0.5) +
  labs(title = "unweighted vs weighted",
       x = "mean expression",
       y = "unweighted rank - weighted rank")

print(rank1)

#check if genes with spatial variation ranked higher compared to genes without spatial variation
rank_unw <- data.frame(
  y = rowData(spe_unweighted)$rank,
  x = rowData(spe_unweighted)$mean,
  top_change_w = top_20_weighted_tf
  
) |> 
  ggplot() +
  geom_point(aes(x = x, y = y, col = top_change_w)) * (blend("lighten") + blend("multiply", alpha = 0.5)) +
  labs(
    x = "mean of logcounts",
    y = "rank",
    title = "unweighted nnSVG"
  )

rank_w <- data.frame(
  y = weighted_rank,
  x = rowData(spe_weighted)$weighted_mean,
  top_change_w = top_20_weighted_tf
) |> 
  ggplot() +
  geom_point(aes(x = x, y = y, col = top_change_w)) * (blend("lighten") + blend("multiply", alpha = 0.5)) +
  labs(
    x = "mean of logcounts",
    y = "rank",
    title = "weighted nnSVG"
  )

wrap_plots(rank_unw, rank_w, nrow=1, guides = "collect") +
  plot_annotation(title = "ranks colored by top 20 changed ranks",
                  caption = "no data cut off")

#check if genes with spatial variation ranked higher compared to genes without spatial variation
#colored by estimated sigma.sq
rank_unw <- data.frame(
  y = rowData(spe_unweighted)$rank,
  x = rowData(spe_unweighted)$mean,
  est_sigma.sq = rowData(spe_unweighted)$sigma.sq
  
) |> 
  ggplot() +
  geom_point(aes(x = x, y = y, col = est_sigma.sq)) +
  scale_color_viridis() +
  labs(
    x = "mean of logcounts",
    y = "rank",
    title = "unweighted nnSVG"
  )

rank_cw <- data.frame(
  y = weighted_rank,
  x = rowData(spe_weighted)$weighted_mean,
  est_sigma.sq = rowData(spe_weighted)$weighted_sigma.sq
) %>% filter(est_sigma.sq < 5) |> 
  ggplot() +
  geom_point(aes(x = x, y = y, col = est_sigma.sq)) +
  scale_color_viridis() +
  labs(
    x = "mean of logcounts",
    y = "rank",
    title = "weighted nnSVG"
  )

wrap_plots(rank_unw, rank_cw, nrow=2) +
  plot_annotation(title = "ranks colored by estimated sigma.sq",
                  caption = "49 genes with est sigma.sq > 5 filtered out from w")

#ridge plots of mean log counts vs rank
#quantiles by means 
#https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009954

#overlay unweighted and constrained weighted
df_unw <- data.frame(
  rank = rowData(spe_unweighted)$rank,
  mean = rowData(spe_unweighted)$mean,
  method = rep("unw", 1853) 
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_cw <- data.frame(
  rank = weighted_rank,
  mean = rowData(spe_weighted)$weighted_mean,
  method = rep("cw", 1853) 
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df <- rbind(df_unw, df_cw) %>% 
  mutate(quantile = as.factor(quantile))

rank_overlay <- ggplot(df, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = method), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - unw & w mean of logcounts",
    x = "rank",
    title = "overlayed unw and constrained w nnSVG"
  )

print(rank_overlay)

#overlay unweighted and constrained weighted using only unw mean
df_unw_2 <- data.frame(
  rank = rowData(spe_unweighted)$rank,
  mean = rowData(spe_unweighted)$mean,
  method = rep("unw", 1853) 
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_cw_2 <- data.frame(
  rank = weighted_rank,
  mean = rowData(spe_unweighted)$mean,
  method = rep("cw", 1853) 
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>%
  tibble::rownames_to_column()

df_2 <- rbind(df_unw_2, df_cw_2) %>% 
  mutate(quantile = as.factor(quantile))

rank_overlay_2 <- ggplot(df_2, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = method), rel_min_height = 0.02, alpha = 0.3) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - unw mean of logcounts",
    x = "rank",
    title = "overlayed unw and constrained w nnSVG"
  )

print(rank_overlay_2)

#ridge plots separated by noise and signal for unw and constrained

frac <- round(dim(spe_unweighted)[1]*0.1*0.1)

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
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - unw mean of logcounts",
    x = "rank",
    title = "unweighted nnSVG"
  )

df_cw_signal <- df_cw %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "signal")

indices <- as.integer(df_cw_signal$rowname)

df_cw_background <- df_cw[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "background")

df <- rbind(df_cw_signal, df_cw_background)

rank_separated_cw <- ggplot(df, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - w mean of logcounts",
    x = "rank",
    title = "constrained weighted nnSVG"
  )

wrap_plots(rank_separated_unw, rank_separated_cw, nrow=1, guides = "collect") +
  plot_annotation(caption = "quantile 1 = lowest means, used top 10% ranks per q")

#ridge plots separated by noise and signal for unw and constrained using unw mean
frac <- round(dim(spe_unweighted)[1]*0.1*0.1)

df_unw_signal_2 <- df_unw_2 %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "signal")

indices <- as.integer(df_unw_signal_2$rowname)

df_unw_background_2 <- df_unw_2[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "background")

df_2 <- rbind(df_unw_signal_2, df_unw_background_2)

rank_separated_unw <- ggplot(df_2, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - unw mean of logcounts",
    x = "rank",
    title = "unweighted nnSVG"
  )

df_cw_signal_2 <- df_cw_2 %>%
  mutate(quantile = as.factor(quantile)) %>%
  group_by(quantile) %>%
  slice_min(order_by = rank, n = frac) %>%
  mutate(grp = "signal")

indices <- as.integer(df_cw_signal_2$rowname)

df_cw_background_2 <- df_cw_2[-indices,] %>%
  mutate(quantile = as.factor(quantile)) %>%
  mutate(grp = "background")

df_2 <- rbind(df_cw_signal_2, df_cw_background_2)

rank_separated_cw <- ggplot(df_2, aes(x = rank, y = quantile)) +
  geom_density_ridges2(aes(fill = grp), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - unw mean of logcounts",
    x = "rank",
    title = "constrained weighted nnSVG"
  )

wrap_plots(rank_separated_unw, rank_separated_cw, nrow=1, guides = "collect") +
  plot_annotation(caption = "quantile 1 = lowest means, used top 10% ranks per q")


#ridge plots for sigma.sq and tau.sq
s_unw <- data.frame(
  sigma.sq = rowData(spe_unweighted)$sigma.sq,
  mean = rowData(spe_unweighted)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = sigma.sq, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "est. sigma^2",
    title = "unweighted nnSVG"
  )

s_cw <- data.frame(
  sigma.sq = rowData(spe_weighted)$weighted_sigma.sq,
  mean = rowData(spe_weighted)$weighted_mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = sigma.sq, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "est. sigma^2",
    title = "con. weighted nnSVG"
  )

wrap_plots(s_unw, s_cw, nrow=1) +
  plot_annotation(title = "sigma.sq ridge plots",
                  caption = "quantile 1 = lowest means, no data filtered out")

t_unw <- data.frame(
  tau.sq = rowData(spe_unweighted)$tau.sq,
  mean = rowData(spe_unweighted)$mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = tau.sq, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "est. tau^2",
    title = "unweighted nnSVG"
  )

t_cw <- data.frame(
  tau.sq = rowData(spe_weighted)$weighted_tau.sq,
  mean = rowData(spe_weighted)$weighted_mean
) %>% mutate(quantile = findInterval(mean, 
                                     quantile(mean, probs=0:9/10))) %>% 
  mutate(quantile = as.factor(quantile)) |> 
  ggplot() +
  geom_density_ridges2(aes(x = tau.sq, y = quantile), rel_min_height = 0.02) +
  theme_ridges(grid = TRUE) +
  labs(
    y = "quantile - mean of logcounts",
    x = "est. tau^2",
    title = "con. weighted nnSVG"
  )

wrap_plots(t_unw, t_cw, nrow=1) +
  plot_annotation(title = "tau.sq ridge plots",
                  caption = "quantile 1 = lowest means, no data filtered out")


#weight estimation figure for HPC

objs <- readRDS(file = "mean_var_project/HPC_V10B01_086_B1_BRISC_estimation_spline.rds")
s_g <- objs[[1]]
r_tilda <- objs[[2]]

dat <- data.frame(
  y = sqrt(s_g),
  x = r_tilda
) 

p1 <- ggplot(dat, aes(x=x, y=y)) +
  geom_point(size=0.5, alpha = 0.2) +
  ylim(c(0, 1.05)) + xlim(c(-1,4)) + 
  labs(
    x = "average log2(count size)",
    y = "sqrt(std dev)"
  ) +
  theme_bw()

#add trend label

p2 <- ggplot(dat, aes(x=x, y=y)) +
  geom_point(size=0.5, alpha = 0.2) +
  geom_smooth(method = 'lm', formula = y ~ splines::bs(x, knots = 10),
              aes(fill = after_scale(color)), alpha = 0.2) +
  ylim(c(0, 1.05)) + xlim(c(-1,4)) + 
  labs(
    x = "average log2(count size)",
    y = "sqrt(std dev)"
  )  +
  annotate("text", x=2.9, y=0.96, label= "gene-wise mean-var trend") +
  theme_bw() 

#add lines?

p3 <- ggplot(dat, aes(x=x, y=y)) +
  geom_smooth(method = 'lm', formula = y ~ splines::bs(x, knots = 10),
              aes(fill = after_scale(color)), alpha = 0.2) +
  ylim(c(0, 1.05)) + xlim(c(-1,4)) + 
  labs(
    x = "fitted log2(count size)",
    y = "sqrt(std dev)"
  )  +
  annotate("text", x=2.9, y=0.96, label= "gene-wise mean-var trend") +
  theme_bw() 

wrap_plots(p1, p2, p3, ncol=3, nrow=1)

mean_data <- data.frame(unweighted = rowData(spe_unweighted)$mean,
                        weighted = rowData(spe_weighted)$weighted_mean,
                        raw = rowMeans(logcounts(spe_unweighted)))


ggpairs(mean_data, aes(alpha = 0.4)) +
  theme(axis.text.x = element_text(angle=90))

dev.off()
