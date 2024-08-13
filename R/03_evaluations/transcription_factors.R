library(SpatialExperiment)
library(dplyr)
library(here)
library(ggplot2)

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


# vis cancer related genes
# Breast cancer related genes
spe_unweighted <- readRDS(here("outputs", "results", "spe_humanBreast_nnSVG.rds"))
genes <- c("ATXN7", "E2F3", "MFHAS1", "SAAL1", "DDX10", "TLN2", "ITCH", "ZHX3")
indices <- which(rowData(spe_unweighted)$gene_name %in% genes)

pdf(here("plots", "transcription_factors_Breast.pdf"))
for (ix in indices) {
  df <- as.data.frame(cbind(spatialCoords(spe_unweighted), expr = logcounts(spe_unweighted)[ix, ]))
  print(ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                       color = expr)) + 
          geom_point(size = 2) + 
          coord_fixed() + 
          scale_y_reverse() + 
          scale_color_viridis_c(name = "logcounts") + 
          ggtitle(rowData(spe_unweighted)$gene_name[ix]) + 
          theme_bw() + 
          theme(plot.title = element_text(face = "italic"), 
                panel.grid = element_blank(), 
                axis.title = element_blank(), 
                axis.text = element_blank(), 
                axis.ticks = element_blank())
  )
}
dev.off()

# Ovarian cancer related genes
spe_unweighted <- readRDS(here("outputs", "results", "spe_humanOvarian_nnSVG.rds"))
genes <- c("TUFT1", "DDX39B")
indices <- which(rowData(spe_unweighted)$gene_name %in% genes)

pdf(here("plots", "transcription_factors_Ovarian.pdf"))
for (ix in indices) {
  df <- as.data.frame(cbind(spatialCoords(spe_unweighted), expr = logcounts(spe_unweighted)[ix, ]))
  print(ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                       color = expr)) + 
          geom_point(size = 3) + 
          coord_fixed() + 
          scale_y_reverse() + 
          scale_color_viridis_c(name = "logcounts") + 
          ggtitle(rowData(spe_unweighted)$gene_name[ix]) + 
          theme_bw() + 
          theme(plot.title = element_text(face = "italic"), 
                panel.grid = element_blank(), 
                axis.title = element_blank(), 
                axis.text = element_blank(), 
                axis.ticks = element_blank())
  )
}
dev.off()

# Lobular Breast cancer related genes
spe_unweighted <- readRDS(here("outputs", "results", "spe_humanLobularBreast_nnSVG.rds"))
genes <- c("PDK1", "PPAT", "ATP6AP1L", "LYRM7", "TRIM35", "TRIM32", "FRAT1", "TIGAR", "ZNF439", "MOSPD2", "BRWD3")
indices <- which(rowData(spe_unweighted)$gene_name %in% genes)

pdf(here("plots", "transcription_factors_LobularBreast.pdf"))
for (ix in indices) {
  df <- as.data.frame(cbind(spatialCoords(spe_unweighted), expr = logcounts(spe_unweighted)[ix, ]))
  print(ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                       color = expr)) + 
          geom_point(size = 2) + 
          coord_fixed() + 
          scale_y_reverse() + 
          scale_color_viridis_c(name = "logcounts") + 
          ggtitle(rowData(spe_unweighted)$gene_name[ix]) + 
          theme_bw() + 
          theme(plot.title = element_text(face = "italic"), 
                panel.grid = element_blank(), 
                axis.title = element_blank(), 
                axis.text = element_blank(), 
                axis.ticks = element_blank())
  )
}
dev.off()

# vis cancer related genes on counts scale
# Breast cancer related genes
spe_unweighted <- readRDS(here("outputs", "results", "spe_humanBreast_nnSVG.rds"))
genes <- c("ATXN7", "E2F3", "MFHAS1", "SAAL1", "DDX10", "TLN2", "ITCH", "ZHX3")
indices <- which(rowData(spe_unweighted)$gene_name %in% genes)

pdf(here("plots", "transcription_factors_Breast_counts.pdf"))
for (ix in indices) {
  df <- as.data.frame(cbind(spatialCoords(spe_unweighted), expr = counts(spe_unweighted)[ix, ]))
  print(ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                       color = expr)) + 
          geom_point(size = 2) + 
          coord_fixed() + 
          scale_y_reverse() + 
          scale_color_viridis_c(name = "counts") + 
          ggtitle(rowData(spe_unweighted)$gene_name[ix]) + 
          theme_bw() + 
          theme(plot.title = element_text(face = "italic"), 
                panel.grid = element_blank(), 
                axis.title = element_blank(), 
                axis.text = element_blank(), 
                axis.ticks = element_blank())
  )
}
dev.off()

# Ovarian cancer related genes
spe_unweighted <- readRDS(here("outputs", "results", "spe_humanOvarian_nnSVG.rds"))
genes <- c("TUFT1", "DDX39B")
indices <- which(rowData(spe_unweighted)$gene_name %in% genes)

pdf(here("plots", "transcription_factors_Ovarian_counts.pdf"))
for (ix in indices) {
  df <- as.data.frame(cbind(spatialCoords(spe_unweighted), expr = counts(spe_unweighted)[ix, ]))
  print(ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                       color = expr)) + 
          geom_point(size = 3) + 
          coord_fixed() + 
          scale_y_reverse() + 
          scale_color_viridis_c(name = "counts") + 
          ggtitle(rowData(spe_unweighted)$gene_name[ix]) + 
          theme_bw() + 
          theme(plot.title = element_text(face = "italic"), 
                panel.grid = element_blank(), 
                axis.title = element_blank(), 
                axis.text = element_blank(), 
                axis.ticks = element_blank())
  )
}
dev.off()

# Lobular Breast cancer related genes
spe_unweighted <- readRDS(here("outputs", "results", "spe_humanLobularBreast_nnSVG.rds"))
genes <- c("PDK1", "PPAT", "ATP6AP1L", "LYRM7", "TRIM35", "TRIM32", "FRAT1", "TIGAR", "ZNF439", "MOSPD2", "BRWD3")
indices <- which(rowData(spe_unweighted)$gene_name %in% genes)

pdf(here("plots", "transcription_factors_LobularBreast_counts.pdf"))
for (ix in indices) {
  df <- as.data.frame(cbind(spatialCoords(spe_unweighted), expr = counts(spe_unweighted)[ix, ]))
  print(ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                       color = expr)) + 
          geom_point(size = 2) + 
          coord_fixed() + 
          scale_y_reverse() + 
          scale_color_viridis_c(name = "counts") + 
          ggtitle(rowData(spe_unweighted)$gene_name[ix]) + 
          theme_bw() + 
          theme(plot.title = element_text(face = "italic"), 
                panel.grid = element_blank(), 
                axis.title = element_blank(), 
                axis.text = element_blank(), 
                axis.ticks = element_blank())
  )
}
dev.off()
