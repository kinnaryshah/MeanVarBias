library(anndata)
library(dplyr)
adata <- read_h5ad("large_scale_sims.h5ad")

library(SpatialExperiment)
counts <- as.matrix(adata$X)
coords <- as.matrix(adata$obs[, c("X_coord", "Y_coord")])
# add one to all the coordinates to make them 1-indexed
coords <- coords + 1

# create SpatialExperiment object
spe <- SpatialExperiment(
  assays = list(counts=t(counts)),
  spatialCoords = coords)

fivenum(rowMeans(counts(spe)))

# run nnSVG

### RUN NNSVG

library(nnSVG)
library(SpatialExperiment)
library(BRISC)
library(BiocParallel)
library(scuttle)
library(Matrix)
library(spoon)

print(dim(spe))

spe <- logNormCounts(spe)

runtime <- system.time({
  tryCatch({
    spe <- nnSVG(
      spe,
      X = NULL,
      assay_name = "logcounts",
      n_neighbors = 10,
      order = "AMMD",
      n_threads = 10,
      verbose = FALSE
    )
  }, error = function(e) {
    print(paste("Error:", e$message))
  })
})

file = paste0("spe_nnSVG_","large",".rds")
saveRDS(spe, file = file)

### run weighted nnSVG

spe <- SpatialExperiment(
  assays = list(counts=t(counts)),
  spatialCoords = coords)

spe <- logNormCounts(spe)

weights <- generate_weights(input = spe,
                            stabilize = TRUE,
                            BPPARAM = MulticoreParam(workers = 10,
                                                     RNGseed = 4))

spe <- weighted_nnSVG(input = spe,
                      w = weights,
                      BPPARAM = MulticoreParam(workers = 10, RNGseed = 5))

file = paste0("spe_weighted_nnSVG_","large",".rds")
saveRDS(spe, file = file)

sum(is.na(rowData(spe)$weighted_LR_stat))

# make ridge plots
library(ggplot2)
library(ggridges)
library(patchwork)

spe_unweighted <- readRDS("spe_nnSVG_large.rds")
spe_weighted <- readRDS("spe_weighted_nnSVG_large.rds")

#ridge plots separated by noise and signal for weighted with unweighted mean
frac <- round(dim(spe_weighted)[1]*0.1*0.1)
n_genes <- dim(spe_weighted)[1]

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
    y = "quaile - unw mean of logcounts",
    x = "rank",
    title = "weighted combined"
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
    y = "decile - unw mean of logcounts",
    x = "rank",
    title = "weighted signal"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw() 

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
    y = "quantile - unw mean of logcounts",
    x = "rank",
    title = "unweighted combined"
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
    y = "decile - unw mean of logcounts",
    x = "rank",
    title = "unweighted signal"
  ) +
  guides(fill=guide_legend(title="group")) +
  coord_cartesian(xlim = c(1, n_genes)) +
  theme_bw()

pdf("visualization_ridge_large.pdf", wid = 21, height = 20)
wrap_plots(list(weighted_ridge_unw_mean, unweighted_ridge), ncol = 2, guides = 'collect')
wrap_plots(list(weighted_ridge_combined_unw_mean, unweighted_ridge_combined), ncol = 2, guides = 'collect')
dev.off()



