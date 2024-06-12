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

# what are the differences between the true negatives in unweighted and weighted
# what genes are we missing with the weighted version? was the weight cut off???
tn_unweighted <- rowData(spe_unweighted)$ground_truth_sigma.sq == 0 & rowData(spe_unweighted)$padj > 0.05
tn_weighted <- rowData(spe_weighted)$ground_truth_sigma.sq == 0 & rowData(spe_weighted)$weighted_padj > 0.05
sum(tn_unweighted, na.rm = TRUE)
sum(tn_weighted, na.rm = TRUE)
# 469, but 9% of observations had weights adjusted

fp_weighted <- rowData(spe_weighted)$ground_truth_sigma.sq == 0 & rowData(spe_weighted)$weighted_padj <= 0.05

# make violin plot of the means of the true negatives and false positives for weighted, and true negatives for unweighted
df_unweighted <- data.frame(mean=rowData(spe_unweighted)$mean[tn_unweighted], method="unweighted")
df_weighted <- data.frame(mean=rowData(spe_unweighted)$mean[tn_weighted], method="weighted")
df_fp <- data.frame(mean=rowData(spe_unweighted)$mean[fp_weighted], method="weighted false positive")
df <- rbind(df_unweighted, df_weighted, df_fp)

ggplot(df, aes(x=method, y=mean)) + geom_violin() + geom_boxplot(width=0.1) + theme_minimal() + labs(title="Mean of true negatives and missed negatives (FDR=0.05)", x="Method", y="Mean")

# read in corresponding weight matrix
weights <- readRDS(paste0(file_dir, "/weights_1.rds"))

# Function to find the modes in order
find_modes <- function(x) {
  unique_x <- unique(x)
  tabulate_x <- tabulate(match(x, unique_x))
  sorted_indices <- order(tabulate_x, decreasing = TRUE)
  sorted_unique_x <- unique_x[sorted_indices]
  sorted_frequencies <- tabulate_x[sorted_indices]
  return(list(values = sorted_unique_x, frequencies = sorted_frequencies))
}

# Find the modes of the values in the weights matrix
# the first and second modes are the constrained weights on the upper and lower end
modes <- find_modes(weights)

# Get the first mode
first_mode <- modes$values[1]
first_mode_count <- modes$frequencies[1]

# Get the second mode, if it exists
if (length(modes$values) > 1) {
  second_mode <- modes$values[2]
  second_mode_count <- modes$frequencies[2]
} else {
  second_mode <- NA
  second_mode_count <- NA
}

first_mode
first_mode_count
second_mode
second_mode_count

# how many genes have weights that match either mode
# genes are on the rows of the weights matrix

first_mode_counts <- numeric(nrow(weights))
second_mode_counts <- numeric(nrow(weights))

# Iterate over each row (gene) in the weights matrix
for (i in 1:nrow(weights)) {
  gene_weights <- weights[i, ]
  first_mode_counts[i] <- sum(gene_weights == first_mode)
}

for (i in 1:nrow(weights)) {
  gene_weights <- weights[i, ]
  second_mode_counts[i] <- sum(gene_weights == second_mode)
}

# Display the mode counts for each gene
first_mode_counts
second_mode_counts

# are the false positives (nonSVGs that were missed by weighted) more likely to have mode 1 or mode 2 weights?
first_mode_counts_fp <- first_mode_counts[fp_weighted]
second_mode_counts_fp <- second_mode_counts[fp_weighted]
both_mode_counts_fp <- first_mode_counts_fp + second_mode_counts_fp

# overlay histograms of the mode counts for the false positives over the mode counts for all genes
# separate histograms for first mode, second mode, and both modes
ggplot(data.frame(counts=first_mode_counts), aes(x=counts)) + geom_histogram(binwidth=1, fill="blue", alpha=0.5) + 
  geom_histogram(data=data.frame(counts=first_mode_counts_fp), aes(x=counts), binwidth=1, fill="red", alpha=0.5) + 
  theme_minimal() + labs(title="Counts of low end constraint per gene in weights matrix", x="Count", y="Frequency", caption = "overlayed false positives") +
  xlim(0, 270)

ggplot(data.frame(counts=second_mode_counts), aes(x=counts)) + geom_histogram(binwidth=1, fill="blue", alpha=0.5) +
  geom_histogram(data=data.frame(counts=second_mode_counts_fp), aes(x=counts), binwidth=1, fill="red", alpha=0.5) + 
  theme_minimal() + labs(title="Counts of high end constraint per gene in weights matrix", x="Count", y="Frequency", caption = "overlayed false positives") +
  xlim(0, 270)

ggplot(data.frame(counts=first_mode_counts+second_mode_counts), aes(x=counts)) + geom_histogram(binwidth=1, fill="blue", alpha=0.5) +
  geom_histogram(data=data.frame(counts=both_mode_counts_fp), aes(x=counts), binwidth=1, fill="red", alpha=0.5) + 
  theme_minimal() + labs(title="Counts of both constraints per gene in weights matrix", x="Count", y="Frequency", caption = "overlayed false positives") +
  xlim(0, 270)
dev.off()


