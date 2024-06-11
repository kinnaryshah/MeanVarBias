library(ggplot2)
library(ggridges)
library(patchwork)
library(SpatialExperiment)
library(dplyr)
library(ggpubr)

file_dir <- "reps_968_50_50per_1000_0.2_to_3_0.5_to_9"
spe_unweighted <- readRDS(paste0(file_dir, "/spe_nnSVG_1.rds"))
spe_weighted <- readRDS(paste0(file_dir, "/spe_weighted_nnSVG_1.rds"))
adj_p_values_weighted <- p.adjust(1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df=2), method = "BH")
rowData(spe_weighted)$weighted_padj <- adj_p_values_weighted

# what are the differences between the true positives in unweighted and weighted
# within the subset of SVGs, which genes have high/low mean
# does nnSVG only get it right for high means? maybe if overall mean is too small, nnSVG throws away
tp_unweighted <- which(rowData(spe_unweighted)$ground_truth_sigma.sq != 0 & rowData(spe_unweighted)$padj <= 0.05)
tp_weighted <- which(rowData(spe_weighted)$ground_truth_sigma.sq != 0 & rowData(spe_weighted)$weighted_padj <= 0.05)
length(tp_unweighted)
length(tp_weighted)

# make violin plot of the means of the true positives for only unweighted, only weighted, and overlaps weighted and unweighted
df_unweighted <- data.frame(mean=rowData(spe_unweighted)$mean[tp_unweighted], method="unweighted")
df_weighted <- data.frame(mean=rowData(spe_unweighted)$mean[tp_weighted], method="weighted")
df <- rbind(df_unweighted, df_weighted)

pdf("mean_lengthscale_vis.pdf")
ggplot(df, aes(x=method, y=mean)) + geom_violin() + geom_boxplot(width=0.1) + theme_minimal() + labs(title="Mean of true positives (FDR=0.05)", x="Method", y="Mean")

# make violin plot of the means of the true positives for only unweighted, only weighted, and overlaps weighted and unweighted
df_unweighted <- data.frame(mean=rowData(spe_unweighted)$ground_truth_sigma.sq[tp_unweighted], method="unweighted")
df_weighted <- data.frame(mean=rowData(spe_weighted)$ground_truth_sigma.sq[tp_weighted], method="weighted")
df <- rbind(df_unweighted, df_weighted)

ggplot(df, aes(x=method, y=mean)) + geom_violin() + geom_boxplot(width=0.1) + theme_minimal() + labs(title="Sigma.sq of true positives (FDR=0.05)", x="Method", y="Sigma.sq")
dev.off()

# what are the differences between the true negatives in unweighted and weighted
# what genes are we missing with the weighted version? was the weight cut off???
tn_unweighted <- rowData(spe_unweighted)$ground_truth_sigma.sq == 0 & rowData(spe_unweighted)$padj > 0.05
tn_weighted <- rowData(spe_weighted)$ground_truth_sigma.sq == 0 & rowData(spe_weighted)$weighted_padj > 0.05
sum(tn_unweighted, na.rm = TRUE)
sum(tn_weighted, na.rm = TRUE)
# 469, but 9% of observations had weights adjusted

fivenum(rowData(spe_unweighted)$mean[tn_unweighted], na.rm = TRUE)
fivenum(rowData(spe_unweighted)$mean[tn_weighted], na.rm = TRUE)

fp_weighted <- which(rowData(spe_weighted)$ground_truth_sigma.sq == 0 & rowData(spe_weighted)$weighted_padj <= 0.05)
fivenum(rowData(spe_unweighted)$mean[fp_weighted], na.rm = TRUE)
