library(ggplot2)
library(SpatialExperiment)
library(patchwork)
library(GGally)
library(dplyr)
library(ggridges)
library(here)
library(ggrepel)

My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  legend.title = element_text(size=20),
  legend.text = element_text(size=20),
  title = element_text(size = 20))

cancer_colors <- c("Low Lengthscale" = "#e1275f", 
                   "Low Mean" = "blue", 
                   "No Cancer" = "gold")

# OVARIAN

spe_nnSVG_Ovarian <- readRDS(here("outputs", "results", "spe_humanOvarian_nnSVG.rds"))
spe_weighted_nnSVG_Ovarian <- readRDS(here("outputs", "results", "spe_humanOvarian_weighted_nnSVG.rds"))

n_genes_Ovarian <- dim(spe_weighted_nnSVG_Ovarian)[1]

# remove one instance of "HSPA14" for duplicate row names
idx <- which(rowData(spe_weighted_nnSVG_Ovarian)$gene_name == "HSPA14")[2]

spe_nnSVG_Ovarian <- spe_nnSVG_Ovarian[-idx,]
spe_weighted_nnSVG_Ovarian <- spe_weighted_nnSVG_Ovarian[-idx,]

df_Ovarian <- data.frame(
  unweighted = rowData(spe_nnSVG_Ovarian)$rank,
  weighted = rowData(spe_weighted_nnSVG_Ovarian)$weighted_rank)

rownames(df_Ovarian) <- rowData(spe_weighted_nnSVG_Ovarian)$gene_name

# cancer genes:
# low lengthscale "LSAMP", "DTX3L", "CASC15", "SLC25A6", "CDK9", "MAP3K8", "GCH1", "PCYT2",
# low mean "TUFT1", "DDX39B"

# create a column in the df called "cancer"
# should be true if in the list above and false if not

cancer_genes_ls <- c("LSAMP", "DTX3L", "CASC15", "SLC25A6", "CDK9", "MAP3K8", "GCH1", "PCYT2")
cancer_genes_mean <- c("TUFT1", "DDX39B")

df_Ovarian$cancer <- ifelse(rownames(df_Ovarian) %in% cancer_genes_ls, "Low Lengthscale",
                            ifelse(rownames(df_Ovarian) %in% cancer_genes_mean, "Low Mean", "No Cancer"))

# plot ranks on either axis
Ovarian_ranks <- ggplot(df_Ovarian, aes(x = unweighted, y = weighted, color = cancer)) +
  geom_point(size = 0.3, alpha=0.3) +
  geom_point(data = filter(df_Ovarian, cancer != "No Cancer"), aes(x = unweighted, y = weighted, color = cancer), size = 3) +
  geom_text_repel(data = filter(df_Ovarian, cancer != "No Cancer"), 
                  aes(x = unweighted, y = weighted, label = rownames(filter(df_Ovarian, cancer != "No Cancer"))),
                  size = 6, max.overlaps = 20, color = "black", force = 4) + 
  xlim(c(1, n_genes_Ovarian)) +
  ylim(c(1, n_genes_Ovarian)) +
  geom_abline(slope=1, intercept=0) +
  ggtitle("Ovarian") +
  xlab("Unweighted Rank") +
  ylab("Weighted Rank") +
  labs(color = "Cancer") +
  scale_color_manual(values = cancer_colors) +
  theme_bw() +
  My_Theme

# BREAST

spe_nnSVG_Breast <- readRDS(here("outputs", "results", "spe_humanBreast_nnSVG.rds"))
spe_weighted_nnSVG_Breast <- readRDS(here("outputs", "results", "spe_humanBreast_weighted_nnSVG.rds"))

n_genes_Breast <- dim(spe_weighted_nnSVG_Breast)[1]

# remove one instance of "HSPA14" for duplicate row names
idx <- which(rowData(spe_weighted_nnSVG_Breast)$gene_name == "HSPA14")[2]

spe_nnSVG_Breast <- spe_nnSVG_Breast[-idx,]
spe_weighted_nnSVG_Breast <- spe_weighted_nnSVG_Breast[-idx,]

df_Breast <- data.frame(
  unweighted = rowData(spe_nnSVG_Breast)$rank,
  weighted = rowData(spe_weighted_nnSVG_Breast)$weighted_rank)

rownames(df_Breast) <- rowData(spe_weighted_nnSVG_Breast)$gene_name

# cancer genes:
# low lengthscale "DDX20", "ASXL2", "MGAT5", "ICOS", "TNFRSF21", "NUP43", "CDK13", "NOL6", "DAPK1", "TRIM21",
# "FAM111B", "SPDYC", "ZBTB1", "PEAK1", "RCCD1", "FTO", "ULK2", "TRIM65", "DNMT3B"
# low mean "ATXN7", "E2F3", "MFHAS1", "SAAL1", "DDX10", "TLN2", "ITCH", "ZHX3"

# create a column in the df called "cancer"
# should be true if in the list above and false if not

cancer_genes_ls <- c("DDX20", "ASXL2", "MGAT5", "ICOS", "TNFRSF21", "NUP43", "CDK13", "NOL6", "DAPK1", "TRIM21",
                  "FAM111B", "SPDYC", "ZBTB1", "PEAK1", "RCCD1", "FTO", "ULK2", "TRIM65", "DNMT3B")
cancer_genes_mean <- c("ATXN7", "E2F3", "MFHAS1", "SAAL1", "DDX10", "TLN2", "ITCH", "ZHX3")

df_Breast$cancer <- ifelse(rownames(df_Breast) %in% cancer_genes_ls, "Low Lengthscale",
                            ifelse(rownames(df_Breast) %in% cancer_genes_mean, "Low Mean", "No Cancer"))

# plot ranks on either axis
Breast_ranks <- ggplot(df_Breast, aes(x = unweighted, y = weighted, color = cancer)) +
  geom_point(size = 0.3, alpha=0.3) +
  geom_point(data = filter(df_Breast, cancer != "No Cancer"), aes(x = unweighted, y = weighted, color = cancer), size = 3) +
  geom_text_repel(data = filter(df_Breast, cancer != "No Cancer"), 
                  aes(x = unweighted, y = weighted, label = rownames(filter(df_Breast, cancer != "No Cancer"))),
                  size = 6, max.overlaps = 20, color = "black", force = 4) + 
  xlim(c(1, n_genes_Breast)) +
  ylim(c(1, n_genes_Breast)) +
  geom_abline(slope=1, intercept=0) +
  ggtitle("Ductal Breast") +
  xlab("Unweighted Rank") +
  ylab("Weighted Rank") +
  labs(color = "Cancer") +
  scale_color_manual(values = cancer_colors) +
  theme_bw() +
  My_Theme

# SUBTYPE BREAST

spe_nnSVG_SubtypeBreast <- readRDS(here("outputs", "results", "spe_humanSubtypeBreast_nnSVG.rds"))
spe_weighted_nnSVG_SubtypeBreast <- readRDS(here("outputs", "results", "spe_humanSubtypeBreast_weighted_nnSVG.rds"))

n_genes_SubtypeBreast <- dim(spe_weighted_nnSVG_SubtypeBreast)[1]

# remove one instance of "HSPA14" for duplicate row names
idx1 <- which(rowData(spe_weighted_nnSVG_SubtypeBreast)$gene_name == "HSPA14")[2]
idx2 <- which(rowData(spe_weighted_nnSVG_SubtypeBreast)$gene_name == "CYB561D2")[1]

spe_nnSVG_SubtypeBreast <- spe_nnSVG_SubtypeBreast[-c(idx1,idx2),]
spe_weighted_nnSVG_SubtypeBreast <- spe_weighted_nnSVG_SubtypeBreast[-c(idx1,idx2),]

df_SubtypeBreast <- data.frame(
  unweighted = rowData(spe_nnSVG_SubtypeBreast)$rank,
  weighted = rowData(spe_weighted_nnSVG_SubtypeBreast)$weighted_rank)

rownames(df_SubtypeBreast) <- rowData(spe_weighted_nnSVG_SubtypeBreast)$gene_name

# cancer genes:
# low lengthscale "CASP8", "ATP11B", "E2F3", "EHMT2", "TREM1", "SMC2", "CD82", "BIRC2", "GABARAPL1",
# "CDYL2", "MSI2", "CTDP1", "DDA1", "JAM2", "TAF1", "FGF13"
# low mean NULL

# create a column in the df called "cancer"
# should be true if in the list above and false if not

cancer_genes_ls <- c("CASP8", "ATP11B", "E2F3", "EHMT2", "TREM1", "SMC2", "CD82", "BIRC2", "GABARAPL1",
                     "CDYL2", "MSI2", "CTDP1", "DDA1", "JAM2", "TAF1", "FGF13")

df_SubtypeBreast$cancer <- ifelse(rownames(df_SubtypeBreast) %in% cancer_genes_ls, "Low Lengthscale", "No Cancer")

# plot ranks on either axis
SubtypeBreast_ranks <- ggplot(df_SubtypeBreast, aes(x = unweighted, y = weighted, color = cancer)) +
  geom_point(size = 0.3, alpha=0.3) +
  geom_point(data = filter(df_SubtypeBreast, cancer != "No Cancer"), aes(x = unweighted, y = weighted, color = cancer), size = 3) +
  geom_text_repel(data = filter(df_SubtypeBreast, cancer != "No Cancer"), 
                  aes(x = unweighted, y = weighted, label = rownames(filter(df_SubtypeBreast, cancer != "No Cancer"))),
                  size = 6, max.overlaps = 20, color = "black", force = 4) + 
  xlim(c(1, n_genes_SubtypeBreast)) +
  ylim(c(1, n_genes_SubtypeBreast)) +
  geom_abline(slope=1, intercept=0) +
  ggtitle("ER+ Breast") +
  xlab("Unweighted Rank") +
  ylab("Weighted Rank") +
  labs(color = "Cancer") +
  theme_bw() +
  scale_color_manual(values = cancer_colors) +
  theme(legend.position="none") +
  My_Theme

# LOBULAR BREAST

spe_nnSVG_LobularBreast <- readRDS(here("outputs", "results", "spe_humanLobularBreast_nnSVG.rds"))
spe_weighted_nnSVG_LobularBreast <- readRDS(here("outputs", "results", "spe_humanLobularBreast_weighted_nnSVG.rds"))

n_genes_LobularBreast <- dim(spe_weighted_nnSVG_LobularBreast)[1]

# remove one instance of "HSPA14" for duplicate row names
idx1 <- which(rowData(spe_weighted_nnSVG_LobularBreast)$gene_name == "HSPA14")[2]
idx2 <- which(rowData(spe_weighted_nnSVG_LobularBreast)$gene_name == "CYB561D2")[1]

spe_nnSVG_LobularBreast <- spe_nnSVG_LobularBreast[-c(idx1,idx2),]
spe_weighted_nnSVG_LobularBreast <- spe_weighted_nnSVG_LobularBreast[-c(idx1,idx2),]

df_LobularBreast <- data.frame(
  unweighted = rowData(spe_nnSVG_LobularBreast)$rank,
  weighted = rowData(spe_weighted_nnSVG_LobularBreast)$weighted_rank)

rownames(df_LobularBreast) <- rowData(spe_weighted_nnSVG_LobularBreast)$gene_name

# cancer genes:
# low lengthscale "ZNF436", "MTF2", "HAGLR", "SATB1", "CD200", "TRIM59", "UCHL1", "LIFR", "FER",
# "NFIL3", "ENTPD1", "CLEC2D", "FLT1", "CGRRF1", "TYRO3", "ADAM10", "NLRP1", "HAUS5", "MECP2"
# low mean "PDK1", "PPAT", "ATP6AP1L", "LYRM7", "TRIM35", "TRIM32", "FRAT1", "TIGAR", "ZNF439", "MOSPD2", "BRWD3"

# create a column in the df called "cancer"
# should be true if in the list above and false if not

cancer_genes_ls <- c("ZNF436", "MTF2", "HAGLR", "SATB1", "CD200", "TRIM59", "UCHL1", "LIFR", "FER",
                  "NFIL3", "ENTPD1", "CLEC2D", "FLT1", "CGRRF1", "TYRO3", "ADAM10", "NLRP1", "HAUS5", "MECP2")
cancer_genes_mean <- c("PDK1", "PPAT", "ATP6AP1L", "LYRM7", "TRIM35", "TRIM32", "FRAT1", "TIGAR", "ZNF439", "MOSPD2", "BRWD3")

df_LobularBreast$cancer <- ifelse(rownames(df_LobularBreast) %in% cancer_genes_ls, "Low Lengthscale",
                            ifelse(rownames(df_LobularBreast) %in% cancer_genes_mean, "Low Mean", "No Cancer"))

# plot ranks on either axis
LobularBreast_ranks <- ggplot(df_LobularBreast, aes(x = unweighted, y = weighted, color = cancer)) +
  geom_point(size = 0.3, alpha=0.3) +
  geom_point(data = filter(df_LobularBreast, cancer != "No Cancer"), aes(x = unweighted, y = weighted, color = cancer), size = 3) +
  geom_text_repel(data = filter(df_LobularBreast, cancer != "No Cancer"), 
                  aes(x = unweighted, y = weighted, label = rownames(filter(df_LobularBreast, cancer != "No Cancer"))),
                  size = 6, max.overlaps = 20, color = "black", force = 4) + 
  xlim(c(1, n_genes_LobularBreast)) +
  ylim(c(1, n_genes_LobularBreast)) +
  geom_abline(slope=1, intercept=0) +
  ggtitle("Lobular Breast") +
  xlab("Unweighted Rank") +
  ylab("Weighted Rank") +
  labs(color = "Cancer") +
  scale_color_manual(values = cancer_colors) +
  theme_bw() +
  My_Theme

scenarios <- c(
  "Low Mean, Lobular Breast",
  "Low Mean, Ductal Breast",
  "Low Mean, Ovarian",
  "Low Mean, ER+ Breast",
  "Small Lengthscale, Lobular Breast",
  "Small Lengthscale, Ductal Breast",
  "Small Lengthscale, Ovarian",
  "Small Lengthscale, ER+ Breast"
)

scenarios_row1 <- list(
  c(A = 20, "A&B" = 11), # 2: 31 genes in List 1, 11 overlaps
  c(A = 18, "A&B" = 8),   # 3: 26 genes in List 1, 8 overlaps
  c(A = 5, "A&B" = 2),  # 1: 7 genes in List 1, 2 overlaps
  c(A = 1, "A&B" = 0)      # 4: 1 gene in List 1, 0 overlaps
)

scenarios_row2 <- list(
  c(A = 65, "A&B" = 19),  # 2: 84 genes in List 1, 19 overlaps
  c(A = 59, "A&B" = 19),  # 3: 78 genes in List 1, 19 overlaps
  c(A = 70, "A&B" = 7),   # 1: 77 genes in List 1, 7 overlaps
  c(A = 43, "A&B" = 16)   # 4: 59 genes in List 1, 16 overlaps
)

create_df <- function(scenarios_data, scenario_names) {
  data <- do.call(rbind, lapply(seq_along(scenarios_data), function(i) {
    scenario <- scenarios_data[[i]]
    data.frame(
      Scenario = scenario_names[i],
      Component = names(scenario),
      Count = as.numeric(scenario)
    )
  }))
  return(data)
}

df_row1 <- create_df(scenarios_row1, scenarios[1:4])
df_row2 <- create_df(scenarios_row2, scenarios[5:9])

df_combined$Component <- ifelse(
  df_combined$Scenario %in% scenarios[1:4] & df_combined$Component == "A&B", "A&B_LowMean",
  ifelse(df_combined$Scenario %in% scenarios[5:8] & df_combined$Component == "A&B", "A&B_LowLengthscale", df_combined$Component)
)

p <- ggplot(df_combined, aes(x = Scenario, y = Count, fill = Component)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c(
      "A" = "#FFD700",              
      "A&B_LowMean" = "#1f77b4",   
      "A&B_LowLengthscale" = "#e1275f" 
    ),
    labels = c(
      "A" = "All Genes",
      "A&B_LowMean" = "Intersection (Low Mean)",
      "A&B_LowLengthscale" = "Intersection (Low Lengthscale)"
    )
  ) +
  labs(
    title = "",
    x = "Dataset",
    y = "Count",
    fill = "Classification"
  ) +
  theme_bw() +
  coord_flip() +
  My_Theme +
  theme(axis.text.y = element_text(size = 12))

# Define layout areas
layout <- c(
  area(t = 1, l = 1, b = 1, r = 2),   
  area(t = 2, l = 1, b = 2, r = 1),  
  area(t = 2, l = 2, b = 2, r = 2), 
  area(t = 3, l = 1, b = 3, r = 1),   
  area(t = 3, l = 2, b = 3, r = 2)   
)

ggsave(here("plots", "main", "compare_ranks_cancer.png"),
       wrap_plots(
         p, 
         LobularBreast_ranks, 
         Breast_ranks,
         Ovarian_ranks, 
         SubtypeBreast_ranks,
         guides = "collect",
         axis_titles = "collect"
       ) + plot_layout(design = layout) +
         plot_annotation(tag_levels = 'A') &
         theme(plot.tag = element_text(size = 25)),
       width = 28, height = 25)
