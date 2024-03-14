library(SpatialExperiment)
library(dplyr)

spe_unweighted <- readRDS("spe_nnSVG.rds")
spe_weighted <- readRDS("spe_weighted_nnSVG.rds")

unweighted_rank <- rowData(spe_unweighted)$rank
weighted_rank <- rowData(spe_weighted)$weighted_rank
true_rank <- rowData(spe_unweighted)$ground_truth_rank

# calculate the correlation between the true rank and the estimated rank
correlation_unweighted <- cor(unweighted_rank, true_rank, method = "spearman")
correlation_weighted <- cor(weighted_rank, true_rank, method = "spearman")

correlation_unweighted
# [1] 0.8799553
correlation_weighted
# [1] 0.8351753

# calculate mean of the ranks of the subset of genes with sigma.sq = 0
# first, get the subset of genes with sigma.sq = 0 from rowData(spe_unweighted$ground_truth_sigma.sq)
# then, get the ranks of these genes from rowData(spe_unweighted$rank)

subset_genes <- which(rowData(spe_unweighted)$ground_truth_sigma.sq == 0)
subset_ranks <- unweighted_rank[subset_genes]
mean(subset_ranks)
# [1] 210.5

subset_ranks <- weighted_rank[subset_genes]
mean(subset_ranks)
# [1] 207.1167

subset_ranks <- true_rank[subset_genes]
mean(subset_ranks)
# [1] 210.5

# get the subset of genes with sigma.sq =/= 0 from rowData(spe_unweighted$ground_truth_sigma.sq)
# then, get the ranks of these genes from rowData(spe_unweighted$rank)

subset_genes <- which(rowData(spe_unweighted)$ground_truth_sigma.sq != 0)
subset_ranks <- unweighted_rank[subset_genes]
mean(subset_ranks)
# [1] 60.5

subset_ranks <- weighted_rank[subset_genes]
mean(subset_ranks)
# [1] 65.575

subset_ranks <- true_rank[subset_genes]
mean(subset_ranks)
# [1] 60.5

# define false positives and false negatives
false_positives_unweighted <- sum(unweighted_rank > true_rank)
false_negatives_unweighted <- sum(unweighted_rank < true_rank)

false_positives_weighted <- sum(weighted_rank > true_rank)
false_negatives_weighted <- sum(weighted_rank < true_rank)

false_positives_unweighted; false_positives_weighted
# [1] 153
# [1] 166
false_negatives_unweighted; false_negatives_weighted
# [1] 141
# [1] 131

# define false positives as the number of genes with sigma.sq = 0 that are ranked lower than the 300 - number of genes with sigma.sq = 0

subset_genes <- which(rowData(spe_unweighted)$ground_truth_sigma.sq == 0)
false_positives_unweighted <- sum(unweighted_rank[subset_genes] < 300 - length(subset_genes))
false_positives_weighted <- sum(weighted_rank[subset_genes] < 300 - length(subset_genes))

false_positives_unweighted; false_positives_weighted
# [1] 0
# [1] 9

