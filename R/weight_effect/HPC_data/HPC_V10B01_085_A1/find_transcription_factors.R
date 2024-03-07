library(SpatialExperiment)
library(dplyr)

spe_unweighted <- readRDS("spe_nnSVG.rds")
spe_weighted <- readRDS("spe_weighted_nnSVG.rds")

unweighted_rank <- rowData(spe_unweighted)$rank
weighted_rank <- rowData(spe_weighted)$weighted_rank
unweighted_mean <- rowData(spe_unweighted)$mean

df <- data.frame(
  unweighted_rank = unweighted_rank,
  weighted_rank = weighted_rank,
  mean = unweighted_mean
)

rownames(df) <- rowData(spe_unweighted)$gene_name

# find 25th percentile of "mean" in df
quantile_25 <- quantile(df$mean, probs = 0.25)

# create a column called "low"
# this will be true if "mean" values less than the 25th percentile
df <- df %>% mutate(low = mean < quantile_25)

# find genes that have high "unweighted_rank" and low "weighted_rank" values

# create a column called "high_unweighted_low_weighted"
# this will be true if "unweighted_rank" values are greater than 0.75*n_genes
# and if "weighted_rank" values are less than 0.25*n_genes
df <- df %>% mutate(high_unweighted_low_weighted = unweighted_rank > 0.75*dim(spe_unweighted)[1] & weighted_rank < 0.25*dim(spe_weighted)[1])

# find genes that have true values for both "low" and "high_unweighted_low_weighted"
# these genes are likely to be transcription factors
transcription_factors <- df %>% filter(low & high_unweighted_low_weighted)
