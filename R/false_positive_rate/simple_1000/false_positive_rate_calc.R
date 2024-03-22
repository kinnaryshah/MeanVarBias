library(SpatialExperiment)

# define the false positive rate and see if it is lower in the weighted version

# 1. define false positives based on the p-values and 
# false positive rate as the total false positives / total genes in the simulation

# 2. compare the false positive rate in the weighted and unweighted version

# 3. plot the false positive rate for the weighted and unweighted version

spe_unweighted <- readRDS("spe_nnSVG.rds")
spe_weighted <- readRDS("spe_weighted_nnSVG.rds")

# 1. define false positives based on the p-values and
# false positive rate as the total false positives / total genes in the simulation

# unweighted

pvalues_unweighted <- rowData(spe_unweighted)$pval
false_positives_unweighted <- sum(pvalues_unweighted < 0.05)
total_genes_unweighted <- dim(spe_unweighted)[1]
false_positive_rate_unweighted <- false_positives_unweighted / total_genes_unweighted

# weighted

pvalues_weighted <- 1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df = 2)
false_positives_weighted <- sum(pvalues_weighted < 0.05)
total_genes_weighted <- dim(spe_weighted)[1]
false_positive_rate_weighted <- false_positives_weighted / total_genes_weighted

# 2. compare the false positive rate in the weighted and unweighted version 

false_positive_rate_unweighted # 0.019
false_positive_rate_weighted # 0.055

