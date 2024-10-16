library(ggplot2)
library(ggridges)
library(patchwork)
library(SpatialExperiment)
library(dplyr)
library(ggpubr)
library(ggforce)

file_dir <- "reps_968_50_50per_1000_0.2_to_3_0.5_to_9"
spe_unweighted <- readRDS(paste0(file_dir, "/spe_nnSVG_1.rds"))
spe_weighted <- readRDS(paste0(file_dir, "/spe_weighted_nnSVG_1_new.rds"))
adj_p_values_weighted <- p.adjust(1 - pchisq(rowData(spe_weighted)$weighted_LR_stat, df=2), method = "BH")
rowData(spe_weighted)$weighted_padj <- adj_p_values_weighted

# what are the differences between the true positives in unweighted and weighted
# within the subset of SVGs, which genes have high/low mean
# does nnSVG only get it right for high means? maybe if overall mean is too small, nnSVG throws away
tp_unweighted <- which(rowData(spe_unweighted)$ground_truth_sigma.sq != 0 & rowData(spe_unweighted)$padj <= 0.05)
tp_weighted <- which(rowData(spe_weighted)$ground_truth_sigma.sq != 0 & rowData(spe_weighted)$weighted_padj <= 0.05)
fn_unweighted <- which(rowData(spe_unweighted)$ground_truth_sigma.sq != 0 & rowData(spe_unweighted)$padj > 0.05)
fn_weighted <- which(rowData(spe_weighted)$ground_truth_sigma.sq != 0 & rowData(spe_weighted)$weighted_padj > 0.05)
length(tp_unweighted)
length(tp_weighted)

# make violin plot of the means of the true positives for only unweighted, only weighted TP and FN
df_unweighted <- data.frame(mean=rowData(spe_unweighted)$mean[tp_unweighted], method="unweighted TP")
df_weighted <- data.frame(mean=rowData(spe_unweighted)$mean[tp_weighted], method="weighted TP")
df_fn_unweighted <- data.frame(mean=rowData(spe_unweighted)$mean[fn_unweighted], method="unweighted FN")
df_fn_weighted <- data.frame(mean=rowData(spe_unweighted)$mean[fn_weighted], method="weighted FN")

df <- rbind(df_unweighted, df_weighted, df_fn_unweighted, df_fn_weighted)

pdf("mean_lengthscale_vis.pdf")
ggplot(df, aes(x=method, y=mean)) + geom_violin() + geom_boxplot(width=0.1) + theme_minimal() + 
  labs(title="Empirical Mean of true positives (FDR=0.05)", x="Method", y="Mean")

# make violin plot of the ground_truth_beta of the true positives for only unweighted, only weighted TP and FN
df_unweighted <- data.frame(mean=rowData(spe_unweighted)$ground_truth_beta[tp_unweighted], method="unweighted TP")
df_weighted <- data.frame(mean=rowData(spe_unweighted)$ground_truth_beta[tp_weighted], method="weighted TP")
df_fn_unweighted <- data.frame(mean=rowData(spe_unweighted)$ground_truth_beta[fn_unweighted], method="unweighted FN")
df_fn_weighted <- data.frame(mean=rowData(spe_unweighted)$ground_truth_beta[fn_weighted], method="weighted FN")

df <- rbind(df_unweighted, df_weighted, df_fn_unweighted, df_fn_weighted)

ggplot(df, aes(x=method, y=mean)) + geom_violin() + geom_boxplot(width=0.1) + theme_minimal() + 
  labs(title="True Beta of true positives (FDR=0.05)", x="Method", y="Beta")

# make violin plot of the sigma.sq of the true positives for only unweighted, only weighted TP and FN
df_unweighted <- data.frame(mean=rowData(spe_unweighted)$ground_truth_sigma.sq[tp_unweighted], method="unweighted TP")
df_weighted <- data.frame(mean=rowData(spe_unweighted)$ground_truth_sigma.sq[tp_weighted], method="weighted TP")
df_fn_unweighted <- data.frame(mean=rowData(spe_unweighted)$ground_truth_sigma.sq[fn_unweighted], method="unweighted FN")
df_fn_weighted <- data.frame(mean=rowData(spe_unweighted)$ground_truth_sigma.sq[fn_weighted], method="weighted FN")

df <- rbind(df_unweighted, df_weighted, df_fn_unweighted, df_fn_weighted)

ggplot(df, aes(x=method, y=mean)) + geom_violin() + geom_boxplot(width=0.1) + theme_minimal() + 
  labs(title="True Sigma.sq of true positives (FDR=0.05)", x="Method", y="Sigma.sq")

# make violin plot of the sigma.sq of the true positives for only unweighted, only weighted TP and FN
df_unweighted <- data.frame(mean=rowData(spe_unweighted)$sigma.sq[tp_unweighted], method="unweighted TP")
df_weighted <- data.frame(mean=rowData(spe_weighted)$weighted_sigma.sq[tp_weighted], method="weighted TP")
df_fn_unweighted <- data.frame(mean=rowData(spe_unweighted)$sigma.sq[fn_unweighted], method="unweighted FN")
df_fn_weighted <- data.frame(mean=rowData(spe_weighted)$weighted_sigma.sq[fn_weighted], method="weighted FN")

df <- rbind(df_unweighted, df_weighted, df_fn_unweighted, df_fn_weighted)

ggplot(df, aes(x=method, y=mean)) + geom_violin() + geom_boxplot(width=0.1) + theme_minimal() + 
  labs(title="Estimated Sigma.sq of true positives (FDR=0.05)", x="Method", y="Sigma.sq")

# make violin plot of the sigma.sq of the true positives for only unweighted TP and FN
df <- rbind(df_unweighted, df_fn_unweighted)

ggplot(df, aes(x=method, y=mean)) + geom_violin() + geom_boxplot(width=0.1) + theme_minimal() + 
  labs(title="Estimated Sigma.sq of true positives (FDR=0.05)", x="Method", y="Sigma.sq")


# make a scatterplot of empirical mean and true sigma.sq SVGs, colored by true positive unweighted, true positive weighted, and neither
# create df of rowData(spe_unweighted)
# add column for true positive unweighted, true positive weighted, and neither
# remove nonSVGs which have ground_truth_sigma.sq = 0
# plot

df <- data.frame(mean=rowData(spe_unweighted)$mean,
                 sigma.sq=rowData(spe_unweighted)$sigma.sq,
                 ground_truth_sigma.sq=rowData(spe_unweighted)$ground_truth_sigma.sq,
                 method="none")

df$method[tp_weighted] <- "weighted"
df$method[tp_unweighted] <- "unweighted"
df$method[intersect(tp_weighted, tp_unweighted)] <- "both"

df <- df[df$ground_truth_sigma.sq != 0,]

ggplot(df, aes(x=mean, y=sigma.sq, color=method)) + geom_point() + theme_minimal() + 
  labs(title="Empirical Mean vs. Estimated Sigma.sq for True SVGs", x="Mean", y="Sigma.sq")

df <- data.frame(mean=rowData(spe_unweighted)$mean,
                 sigma.sq=rowData(spe_unweighted)$sigma.sq,
                 ground_truth_sigma.sq=rowData(spe_unweighted)$ground_truth_sigma.sq,
                 method="none")

df$method[tp_weighted] <- "weighted"
df$method[tp_unweighted] <- "unweighted"
df$method[intersect(tp_weighted, tp_unweighted)] <- "both"

df <- df[df$ground_truth_sigma.sq != 0,]

ggplot(df, aes(x=mean, y=ground_truth_sigma.sq, color=method)) + geom_point() + theme_minimal() + 
  labs(title="Empirical Mean vs. True Sigma.sq for True SVGs", x="Mean", y="Sigma.sq")


# make a scatterplot of the empirical mean and estimated lengthscale, colored by true positive unweighted

df <- data.frame(mean=rowData(spe_unweighted)$mean,
                 lengthscale=log(1/rowData(spe_unweighted)$phi),
                 ground_truth_sigma.sq=rowData(spe_unweighted)$ground_truth_sigma.sq,
                 method="unweighted FN")

df$method[tp_unweighted] <- "unweighted TP"

df <- df[df$ground_truth_sigma.sq != 0,]

ggplot(df, aes(x=mean, y=lengthscale, color=method)) + geom_point() + theme_minimal() + 
  labs(title="Empirical Mean vs. Estimated Unweighted Lengthscale for True SVGs", x="Mean", y="Log Lengthscale", caption = "True lengthscale=50")

# make a scatterplot of the empirical mean and estimated lengthscale, colored by true positive weighted

df <- data.frame(mean=rowData(spe_unweighted)$mean,
                 lengthscale=log(1/rowData(spe_weighted)$weighted_phi),
                 ground_truth_sigma.sq=rowData(spe_weighted)$ground_truth_sigma.sq,
                 method="weighted FN")

df$method[tp_weighted] <- "weighted TP"

df <- df[df$ground_truth_sigma.sq != 0,]

ggplot(df, aes(x=mean, y=lengthscale, color=method)) + geom_point() + theme_minimal() + 
  labs(title="Empirical Mean vs. Estimated Weighted Lengthscale for True SVGs", x="Mean", y="Log Lengthscale", caption = "True lengthscale=50")

# make a scatterplot of the ground_truth_sigma.sq and estimated lengthscale, colored by true positive unweighted

df <- data.frame(ground_truth_sigma.sq=rowData(spe_unweighted)$ground_truth_sigma.sq,
                 lengthscale=log(1/rowData(spe_unweighted)$phi),
                 method="unweighted FN")

df$method[tp_unweighted] <- "unweighted TP"

df <- df[df$ground_truth_sigma.sq != 0,]

ggplot(df, aes(x=ground_truth_sigma.sq, y=lengthscale, color=method)) + geom_point() + theme_minimal() + 
  labs(title="True Sigma.sq vs. Estimated Unweighted Lengthscale for True SVGs", x="Sigma.sq", y="Log Lengthscale", caption = "True lengthscale=50")

# make a scatterplot of the ground_truth_sigma.sq and estimated lengthscale, colored by true positive weighted

df <- data.frame(ground_truth_sigma.sq=rowData(spe_weighted)$ground_truth_sigma.sq,
                 lengthscale=log(1/rowData(spe_weighted)$weighted_phi),
                 method="weighted FN")

df$method[tp_weighted] <- "weighted TP"

df <- df[df$ground_truth_sigma.sq != 0,]

ggplot(df, aes(x=ground_truth_sigma.sq, y=lengthscale, color=method)) + geom_point() + theme_minimal() + 
  labs(title="True Sigma.sq vs. Estimated Weighted Lengthscale for True SVGs", x="Sigma.sq", y="Log Lengthscale", caption = "True lengthscale=50")

# make a scatterplot of the ground_truth_sigma.sq and estimated sigma.sq, colored by true positive unweighted

df <- data.frame(ground_truth_sigma.sq=rowData(spe_unweighted)$ground_truth_sigma.sq,
                 estimated_sigma.sq=rowData(spe_unweighted)$sigma.sq,
                 method="unweighted FN")

df$method[tp_unweighted] <- "unweighted TP"

df <- df[df$ground_truth_sigma.sq != 0,]

ggplot(df, aes(x=ground_truth_sigma.sq, y=estimated_sigma.sq, color=method)) + geom_point() + theme_minimal() + 
  labs(title="True Sigma.sq vs. Estimated Unweighted Sigma.sq for True SVGs", x="True Sigma.sq", y="Estimated Sigma.sq")

# make a scatterplot of the ground_truth_sigma.sq and estimated sigma.sq, colored by true positive weighted

df <- data.frame(ground_truth_sigma.sq=rowData(spe_weighted)$ground_truth_sigma.sq,
                 estimated_sigma.sq=log(rowData(spe_weighted)$weighted_sigma.sq),
                 method="weighted FN")

df$method[tp_weighted] <- "weighted TP"

df <- df[df$ground_truth_sigma.sq != 0,]

ggplot(df, aes(x=ground_truth_sigma.sq, y=estimated_sigma.sq, color=method)) + geom_point() + theme_minimal() + 
  labs(title="True Sigma.sq vs. Log Estimated Weighted Sigma.sq for True SVGs", x="True Sigma.sq", y="Log Estimated Sigma.sq")

# make a scatterplot of the estimated sigma.sq and estimated lengthscale, colored by true positive unweighted

df <- data.frame(sigma.sq=rowData(spe_unweighted)$sigma.sq,
                 lengthscale=log(1/rowData(spe_unweighted)$phi),
                 method="unweighted FN")

df$method[tp_unweighted] <- "unweighted TP"

df <- df[df$sigma.sq != 0,]

ggplot(df, aes(x=sigma.sq, y=lengthscale, color=method)) + geom_point() + theme_minimal() + 
  labs(title="Estimated Unweighted Sigma.sq vs. Estimated Unweighted Lengthscale for True SVGs", x="Estimated Sigma.sq", y="Log Lengthscale", caption = "True lengthscale=50")

# make a scatterplot of the estimated sigma.sq and estimated lengthscale, colored by true positive weighted

df <- data.frame(sigma.sq=log(rowData(spe_weighted)$weighted_sigma.sq),
                 lengthscale=log(1/rowData(spe_weighted)$weighted_phi),
                 method="weighted FN")

df$method[tp_weighted] <- "weighted TP"

df <- df[df$sigma.sq != 0,]

ggplot(df, aes(x=sigma.sq, y=lengthscale, color=method)) + geom_point() + theme_minimal() + 
  labs(title="Log Estimated Weighted Sigma.sq vs. Estimated Weighted Lengthscale for True SVGs", x="Log Estimated Sigma.sq", y="Log Lengthscale", caption = "True lengthscale=50")

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

# make a scatterplot of the empirical mean and estimated lengthscale colored by true negative unweighted

df <- data.frame(mean=rowData(spe_unweighted)$mean[tn_unweighted],
                 lengthscale=log(1/rowData(spe_unweighted)$phi[tn_unweighted]),
                 method="unweighted TN")

ggplot(df, aes(x=mean, y=lengthscale, color=method)) + geom_point() + theme_minimal() +
  labs(title="Empirical Mean vs. Estimated Unweighted Lengthscale for True nonSVGs", x="Empirical Mean", y="Log Lengthscale", caption = "True lengthscale=50")

# make a scatterplot of the empirical mean and estimated lengthscale colored by true negative weighted

df <- data.frame(mean=rowData(spe_unweighted)$mean,
                 lengthscale=log(1/rowData(spe_weighted)$weighted_phi),
                 ground_truth_sigma.sq=rowData(spe_weighted)$ground_truth_sigma.sq,
                 method="weighted FP")

df$method[tn_weighted] <- "weighted TN"

df <- df[df$ground_truth_sigma.sq == 0,]

ggplot(df, aes(x=mean, y=lengthscale, color=method)) + geom_point() + theme_minimal() +
  labs(title="Empirical Mean vs. Estimated Weighted Lengthscale for True nonSVGs", x="Empirical Mean", y="Log Lengthscale", caption = "True lengthscale=50")

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
# genes are on the cols of the weights matrix

first_mode_counts <- numeric(ncol(weights))
second_mode_counts <- numeric(ncol(weights))

# Iterate over each row (gene) in the weights matrix
for (i in 1:ncol(weights)) {
  gene_weights <- weights[, i]
  first_mode_counts[i] <- sum(gene_weights == first_mode)
}

for (i in 1:ncol(weights)) {
  gene_weights <- weights[, i]
  second_mode_counts[i] <- sum(gene_weights == second_mode)
}

# are the false positives (nonSVGs that were missed by weighted) more likely to have mode 1 or mode 2 weights?
first_mode_counts_fp <- first_mode_counts[fp_weighted]
second_mode_counts_fp <- second_mode_counts[fp_weighted]
both_mode_counts_fp <- first_mode_counts_fp + second_mode_counts_fp

# overlay histograms of the mode counts for the false positives over the mode counts for all genes
# separate histograms for first mode, second mode, and both modes
ggplot(data.frame(counts=first_mode_counts), aes(x=counts)) + geom_histogram(binwidth=1, fill="blue", alpha=0.5) + 
  geom_histogram(data=data.frame(counts=first_mode_counts_fp), aes(x=counts), binwidth=1, fill="red", alpha=0.5) + 
  theme_minimal() + labs(title="Counts of low end constraint per gene in weights matrix", x="Count", y="Frequency", caption = "overlayed false positives") +
  facet_zoom(ylim=c(0, 35), xlim=c(0, 600))

ggplot(data.frame(counts=second_mode_counts), aes(x=counts)) + geom_histogram(binwidth=1, fill="blue", alpha=0.5) +
  geom_histogram(data=data.frame(counts=second_mode_counts_fp), aes(x=counts), binwidth=1, fill="red", alpha=0.5) + 
  theme_minimal() + labs(title="Counts of high end constraint per gene in weights matrix", x="Count", y="Frequency", caption = "overlayed false positives") +
  facet_zoom(ylim=c(0, 35), xlim=c(0, 600))

ggplot(data.frame(counts=first_mode_counts+second_mode_counts), aes(x=counts)) + geom_histogram(binwidth=1, fill="blue", alpha=0.5) +
  geom_histogram(data=data.frame(counts=both_mode_counts_fp), aes(x=counts), binwidth=1, fill="red", alpha=0.5) + 
  theme_minimal() + labs(title="Counts of both constraints per gene in weights matrix", x="Count", y="Frequency", caption = "overlayed false positives") +
  facet_zoom(ylim=c(0, 35), xlim=c(0, 600))

dev.off()

# are there any false positives that have no weights cut off
# i.e. all weights are neither of the modes
no_mode_fp <- which(first_mode_counts_fp + second_mode_counts_fp == 0)

# see if cut off weights impact large lengthscale values for the false positives
# x axis is first_mode_counts_fp, y axis is estimated log lengthscale
# colored by whether the gene was a true negative or false positive
# only weighted because the unweighted weights are all 1

df <- data.frame(lengthscale=log(1/rowData(spe_weighted)$weighted_phi),
                 first_mode_counts=first_mode_counts,
                 ground_truth_sigma.sq=rowData(spe_weighted)$ground_truth_sigma.sq,
                 method="weighted FP")
df$method[tn_weighted] <- "weighted TN"

df <- df[df$ground_truth_sigma.sq == 0,]

ggplot(df, aes(x=first_mode_counts, y=lengthscale, color=method)) + geom_point() + theme_minimal() +
  labs(title="Low end constraint vs. Estimated Weighted Lengthscale for True nonSVGs", x="Count Low End Constraint", y="Log Lengthscale", caption = "True lengthscale=50")

# high end constraint

df <- data.frame(lengthscale=log(1/rowData(spe_weighted)$weighted_phi),
                 second_mode_counts=second_mode_counts,
                 ground_truth_sigma.sq=rowData(spe_weighted)$ground_truth_sigma.sq,
                 method="weighted FP")
df$method[tn_weighted] <- "weighted TN"

df <- df[df$ground_truth_sigma.sq == 0,]

ggplot(df, aes(x=second_mode_counts, y=lengthscale, color=method)) + geom_point() + theme_minimal() +
  labs(title="High end constraint vs. Estimated Weighted Lengthscale", x="Count High End Constraint", y="Log Lengthscale", caption = "True lengthscale=50")

# both constraints

df <- data.frame(lengthscale=log(1/rowData(spe_weighted)$weighted_phi),
                 both_mode_counts=first_mode_counts+second_mode_counts,
                 ground_truth_sigma.sq=rowData(spe_weighted)$ground_truth_sigma.sq,
                 method="weighted FP")
df$method[tn_weighted] <- "weighted TN"

df <- df[df$ground_truth_sigma.sq == 0,]

ggplot(df, aes(x=both_mode_counts, y=lengthscale, color=method)) + geom_point() + theme_minimal() +
  labs(title="Both Constraints vs. Estimated Weighted Lengthscale", x="Count Both Constraints", y="Log Lengthscale", caption = "True lengthscale=50")

