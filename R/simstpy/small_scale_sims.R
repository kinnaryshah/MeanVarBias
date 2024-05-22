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
      n_threads = 1,
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
                            BPPARAM = MulticoreParam(workers = 1,
                                                     RNGseed = 4))

spe <- weighted_nnSVG(input = spe,
                      w = weights,
                      BPPARAM = MulticoreParam(workers = 1, RNGseed = 5))

# find rows with all zeroes in matrix called "weighted_logcounts"
zero_rows <- colSums(weighted_logcounts) == 0

file = paste0("spe_weighted_nnSVG_","large",".rds")
saveRDS(spe, file = file)

sum(is.na(rowData(spe)$weighted_LR_stat))

# make ridge plots
library(ggplot2)
library(ggridges)
library(patchwork)

spe_unweighted <- readRDS("spe_nnSVG_1.rds")
spe_weighted <- readRDS("spe_weighted_nnSVG_1.rds")

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
    y = "quantile - unw mean of logcounts",
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

pdf("visualization_ridge_1.pdf", width = 21, height = 20)
wrap_plots(list(weighted_ridge_unw_mean, unweighted_ridge), ncol = 2, guides = 'collect')
wrap_plots(list(weighted_ridge_combined_unw_mean, unweighted_ridge_combined), ncol = 2, guides = 'collect')
dev.off()

# make FDR plots

# create variable called "ground_truth_SVG" that is FALSE if the gene name contains nonSVG and TRUE otherwise
rowData(spe_unweighted)$ground_truth_SVG <- !grepl("nonSVG", rownames(spe_unweighted))
rowData(spe_weighted)$ground_truth_SVG <- !grepl("nonSVG", rownames(spe_unweighted))

spe_unweighted <- readRDS(paste0("spe_nnSVG_", "large", ".rds"))
p1 <- data.frame(pval = rowData(spe_unweighted)$pval) %>%
  ggplot(aes(x = pval)) +
  geom_histogram(binwidth = 0.01) +
  ggtitle(paste0("unweighted p-value distribution for ", "75% nonSVGs")) +
  # add limits to y-axis
  ylim(0, 300) +
  theme_minimal()

spe_weighted <- readRDS(paste0("spe_weighted_nnSVG_", "large", ".rds"))
p2 <- data.frame(pval = 1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df=2)) %>%
  ggplot(aes(x = pval)) +
  geom_histogram(binwidth = 0.01) +
  ggtitle(paste0("weighted p-value distribution for ", "75% nonSVGs")) +
  # add limits to y-axis
  ylim(0, 300) +
  theme_minimal()

pdf(paste0("p_value_distribution_", "large", ".pdf"))
print(p1)
print(p2)
dev.off()

create_FDR_df <- function(spe_unweighted, spe_weighted) {
  
  alpha_levels <- seq(0.001, 0.1, by = 0.001)  
  
  # UNWEIGHTED
  
  adj_p_values_unweighted <- rowData(spe_unweighted)$padj
  
  fdr_percent_unweighted <- sapply(alpha_levels, function(alpha) {
    num_positive_tests <- sum(adj_p_values_unweighted <= alpha)
    # calculate the number of false positives as the number of non DE genes with adjusted p-value <= alpha
    num_false_positives <- sum(rowData(spe_unweighted)$padj[rowData(spe_unweighted)$ground_truth_SVG == FALSE] <= alpha)
    if (num_positive_tests > 0) {
      fdr <- num_false_positives / num_positive_tests * 100
    } else {
      fdr <- 0
    }
    return(fdr)
  })
  
  
  # WEIGHTED
  
  adj_p_values_weighted <- p.adjust(1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df=2), method = "BH")
  rowData(spe_weighted)$weighted_padj <- adj_p_values_weighted
  
  fdr_percent_weighted <- sapply(alpha_levels, function(alpha) {
    num_positive_tests <- sum(adj_p_values_weighted <= alpha, na.rm = TRUE)
    # calculate the number of false positives as the number of non DE genes with adjusted p-value <= alpha
    num_false_positives <- sum(rowData(spe_weighted)$weighted_padj[rowData(spe_weighted)$ground_truth_SVG == FALSE] <= alpha, na.rm = TRUE)
    if (num_positive_tests > 0) {
      fdr <- num_false_positives / num_positive_tests * 100
    } else {
      fdr <- 0
    }
    return(fdr)
  })
  
  df <- data.frame(fdr_percent_unweighted = fdr_percent_unweighted,
                   fdr_percent_weighted = fdr_percent_weighted,
                   alpha = alpha_levels)
  
  return(df)
}
df <- create_FDR_df(spe_unweighted, spe_weighted)
p_FDR <- df %>%
  ggplot(aes(x = alpha)) +
  geom_line(aes(y = fdr_percent_unweighted, color = "unweighted")) +
  geom_line(aes(y = fdr_percent_weighted, color = "weighted")) +
  ggtitle("FDR vs alpha") +
  theme_bw() +
  theme(legend.position = "top")

create_TNR_TPR_df <- function(spe_unweighted, spe_weighted) {
  
  fdr_levels <- seq(0.001, 0.1, by = 0.001)  
  tpr_levels_unweighted <- numeric(length(fdr_levels))
  tnr_levels_unweighted <- numeric(length(fdr_levels))
  tpr_levels_weighted <- numeric(length(fdr_levels))
  tnr_levels_weighted <- numeric(length(fdr_levels))
  
  # UNWEIGHTED
  
  for (k in 1:length(fdr_levels)) {
    tp <- sum(rowData(spe_unweighted)$ground_truth_SVG != FALSE & rowData(spe_unweighted)$padj <= fdr_levels[k])
    tn <- sum(rowData(spe_unweighted)$ground_truth_SVG == FALSE & rowData(spe_unweighted)$padj > fdr_levels[k])
    fn <- sum(rowData(spe_unweighted)$ground_truth_SVG != FALSE & rowData(spe_unweighted)$padj > fdr_levels[k])
    fp <- sum(rowData(spe_unweighted)$ground_truth_SVG == FALSE & rowData(spe_unweighted)$padj <= fdr_levels[k])
    
    tpr_levels_unweighted[k] <- tp / (tp + fn)
    tnr_levels_unweighted[k] <- tn / (tn + fp)
  }
  
  # WEIGHTED
  
  adj_p_values_weighted <- p.adjust(1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df=2), method = "BH")
  rowData(spe_weighted)$weighted_padj <- adj_p_values_weighted
  
  for (k in 1:length(fdr_levels)) {
    tp <- sum(rowData(spe_weighted)$ground_truth_SVG != FALSE & rowData(spe_weighted)$weighted_padj <= fdr_levels[k], na.rm = TRUE)
    tn <- sum(rowData(spe_weighted)$ground_truth_SVG == FALSE & rowData(spe_weighted)$weighted_padj > fdr_levels[k], na.rm = TRUE)
    fn <- sum(rowData(spe_weighted)$ground_truth_SVG != FALSE & rowData(spe_weighted)$weighted_padj > fdr_levels[k], na.rm = TRUE)
    fp <- sum(rowData(spe_weighted)$ground_truth_SVG == FALSE & rowData(spe_weighted)$weighted_padj <= fdr_levels[k], na.rm = TRUE)
    
    tpr_levels_weighted[k] <- tp / (tp + fn)
    tnr_levels_weighted[k] <- tn / (tn + fp)
  }
  
  df <- data.frame(tpr_levels_unweighted = tpr_levels_unweighted,
                   tnr_levels_unweighted = tnr_levels_unweighted,
                   tpr_levels_weighted = tpr_levels_weighted,
                   tnr_levels_weighted = tnr_levels_weighted,
                   fdr_levels = fdr_levels)
  
  return(df)
}
df <- create_TNR_TPR_df(spe_unweighted, spe_weighted)

# plot average TPR and TNR for different FDR levels
p_TPR <- df %>%
  ggplot(aes(x = fdr_levels)) +
  geom_line(aes(y = tpr_levels_unweighted, color = "unweighted")) +
  geom_line(aes(y = tpr_levels_weighted, color = "weighted")) +
  ggtitle("TPR vs FDR") +
  theme_bw() +
  theme(legend.position = "top")

p_TNR <- df %>%
  ggplot(aes(x = fdr_levels)) +
  geom_line(aes(y = tnr_levels_unweighted, color = "unweighted")) +
  geom_line(aes(y = tnr_levels_weighted, color = "weighted")) +
  ggtitle("TNR vs FDR") +
  theme_bw() +
  theme(legend.position = "top")

# save plots
pdf("FDR_TPR_TNR_large.pdf")
print(p_FDR)
print(p_TPR)
print(p_TNR)
dev.off()

# make ridge plots only for ground truth SVGs
library(ggplot2)
library(ggridges)
library(patchwork)

# only keep rows with ground truth SVGs
# we can tell if a gene is a ground truth SVG by finding if the gene name does NOT contain "nonSVG"
indices <- !grepl("nonSVG", rownames(spe_unweighted))
spe_unweighted <- spe_unweighted[indices,]
spe_weighted <- spe_weighted[indices,]

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
    y = "quantile - unw mean of logcounts",
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

pdf("visualization_ridge_SVGs_large.pdf", width = 21, height = 20)
wrap_plots(list(weighted_ridge_unw_mean, unweighted_ridge), ncol = 2, guides = 'collect')
wrap_plots(list(weighted_ridge_combined_unw_mean, unweighted_ridge_combined), ncol = 2, guides = 'collect')
dev.off()

