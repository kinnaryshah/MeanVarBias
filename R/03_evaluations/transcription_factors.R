library(SpatialExperiment)
library(dplyr)
library(here)

spe_unweighted <- readRDS(here("outputs", "results", paste0("spe_human", "Breast", "_nnSVG.rds")))
spe_weighted <- readRDS(here("outputs", "results", paste0("spe_human", "Breast", "_noise_weighted_nnSVG.rds")))

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
# this will be true if "unweighted_rank" values are greater than 0.9*n_genes
# and if "weighted_rank" values are less than 0.1*n_genes
df <- df %>% mutate(high_unweighted_low_weighted = unweighted_rank > 0.9*dim(spe_unweighted)[1] & weighted_rank < 0.1*dim(spe_weighted)[1])

# find genes that have true values for both "low" and "high_unweighted_low_weighted"
# these genes are likely to be transcription factors
transcription_factors <- df %>% filter(low & high_unweighted_low_weighted)

# print out a list of transcription factors separated by commas with no ""
cat(paste0(rownames(transcription_factors), collapse = ", "))
# TMEM39B, ETAA1, ATXN7, BBS7, MFSD8, ETFDH, E2F3, FIG4, TSPYL4, METTL2B, MFHAS1, MAK16, GKAP1, SAAL1, DDX10, B3GLCT, GNB5, TLN2, VPS33B, ASB7, MOSMO, MED26, ZNF227, ITCH, ZHX3, ZBTB21


# make a loop out of the above code to find transcription factors for all tissues in the list dat=c("Breast","Ovarian", "LobularBreast", "SubtypeBreast)

for (tissue in c("Breast","Ovarian", "LobularBreast", "SubtypeBreast")) {
  spe_unweighted <- readRDS(here("outputs", "results", paste0("spe_human", tissue, "_nnSVG.rds")))
  spe_weighted <- readRDS(here("outputs", "results", paste0("spe_human", tissue, "_noise_weighted_nnSVG.rds")))
  
  unweighted_rank <- rowData(spe_unweighted)$rank
  weighted_rank <- rowData(spe_weighted)$weighted_rank
  unweighted_mean <- rowData(spe_unweighted)$mean
  
  df <- data.frame(
    unweighted_rank = unweighted_rank,
    weighted_rank = weighted_rank,
    mean = unweighted_mean
  )
  
  # fix error non-unique value when setting 'row.names' for any genes
  # find non unique genes
  non_unique_genes <- rowData(spe_unweighted)$gene_name[duplicated(rowData(spe_unweighted)$gene_name)]
  
  # change second index of non unique genes to "gene_name_1" in rowData(spe_unweighted)$gene_name
  for (gene in non_unique_genes) {
    index <- which(rowData(spe_unweighted)$gene_name == gene)
    rowData(spe_unweighted)$gene_name[index[2]] <- paste0(gene, "_1")
  }
  
  rownames(df) <- rowData(spe_unweighted)$gene_name
  
  # find 25th percentile of "mean" in df
  quantile_25 <- quantile(df$mean, probs = 0.25)
  
  # create a column called "low"
  # this will be true if "mean" values less than the 25th percentile
  df <- df %>% mutate(low = mean < quantile_25)
  
  # find genes that have high "unweighted_rank" and low "weighted_rank" values
  
  # create a column called "high_unweighted_low_weighted"
  # this will be true if "unweighted_rank" values are greater than 0.9*n_genes
  # and if "weighted_rank" values are less than 0.1*n_genes
  df <- df %>% mutate(high_unweighted_low_weighted = unweighted_rank > 0.9*dim(spe_unweighted)[1] & weighted_rank < 0.1*dim(spe_weighted)[1])
  
  # find genes that have true values for both "low" and "high_unweighted_low_weighted"
  # these genes are likely to be transcription factors
  transcription_factors <- df %>% filter(low & high_unweighted_low_weighted)
  
  # print out a list of transcription factors separated by commas with no ""
  cat(paste0("Tissue: ", tissue, "\n"))
  cat(paste0(rownames(transcription_factors), collapse = ", "))
  cat("\n\n")
}


# Tissue: Breast
# TMEM39B, ETAA1, ATXN7, BBS7, MFSD8, ETFDH, E2F3, FIG4, TSPYL4, METTL2B, MFHAS1, MAK16, GKAP1, SAAL1, DDX10, B3GLCT, GNB5, TLN2, VPS33B, ASB7, MOSMO, MED26, ZNF227, ITCH, ZHX3, ZBTB21
# 
# Tissue: Ovarian
# TUFT1, EHHADH, DDX39B, NAV2, AC026471.4, SMYD4, HEXIM2
# 
# Tissue: LobularBreast
# TMEM51, PDK1, QTRT2, OSBPL11, TMEM44, CPLX1, PPAT, ATP6AP1L, LYRM7, MSH5, ZBTB24, SHPRH, TRIM35, NCALD, SNX30, TRIM32, FRAT1, RAB11FIP2, TIGAR, INTS13, DLEU2, ATXN1L, SNHG30, SEPTIN4, ZNF439, ZNF181, TMEM191B, ZNF74, MOSPD2, FAAH2, BRWD3
# 
# Tissue: SubtypeBreast
# RELA-DT