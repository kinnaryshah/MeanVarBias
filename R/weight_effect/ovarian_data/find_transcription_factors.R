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

# fix error non-unique value when setting 'row.names': ‘HSPA14’ 
# find index of "HSPA14" in rowData(spe_unweighted)$gene_name
index <- which(rowData(spe_unweighted)$gene_name == "HSPA14")

# change second index of "HSPA14" to "HSPA14_1" in rowData(spe_unweighted)$gene_name
rowData(spe_unweighted)$gene_name[index[2]] <- "HSPA14_1"

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
df <- df %>% mutate(high_unweighted_low_weighted = unweighted_rank > 0.9*dim(spe_unweighted)[1] & weighted_rank < 0.1*dim(spe_weighted)[1])

# find genes that have true values for both "low" and "high_unweighted_low_weighted"
# these genes are likely to be transcription factors
transcription_factors <- df %>% filter(low & high_unweighted_low_weighted)

#HEXIM2, SGF29, ARID3B, ESR1, DDX39B are the transcription factors https://www.tfcheckpoint.org/index.php/index.php/index.php/index.php/index.php/index.php/index.php/index.php/index.php/index.php/index.php/search
