library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(here)
library(SpatialExperiment)
library(ggplot2)

# create a function to get a list of genes from the argument spe object
# we want to get the list of genes that are SVGs for the weighted and unweighted versions of the spe object
# if only_small_ls is TRUE, we want the small lengthscale genes only
# if only_small_ls is FALSE, we want all genes
# return a list of genes for the unweighted and weighted versions of the spe object that can be input in GO analysis

get_gene_lists <- function(spe_unweighted, spe_weighted, only_small_ls = F, low, high) {
  if (only_small_ls) {
    idx <- which(rowData(spe_unweighted)$phi > low & rowData(spe_unweighted)$phi < high)
    # subset to small lengthscale genes in both spe
    spe_unweighted <- spe_unweighted[idx,]
    spe_weighted <- spe_weighted[idx,]
  }
  
  # get the gene list for the unweighted spe using rowData(spe_unweighted)$padj < 0.05
  # gene id is in rowData(spe_unweighted)$gene_id
  gene_list_unweighted <- rowData(spe_unweighted)$gene_id[which(rowData(spe_unweighted)$padj < 0.05)]
  
  # get the gene list for the weighted spe using rowData(spe_weighted)$weighted_padj < 0.05
  gene_list_weighted <- rowData(spe_weighted)$gene_id[which(rowData(spe_weighted)$weighted_padj < 0.05)]
  
  # return the gene list for the unweighted and weighted spe
  return(list(gene_list_unweighted, gene_list_weighted))
  
}

dat=c("LobularBreast")

for (d in dat) {
  spe_unweighted <- readRDS(here("outputs", "results", paste0("spe_human", d, "_nnSVG.rds")))
  spe_weighted <- readRDS(here("outputs", "results", paste0("spe_human", d, "_weighted_nnSVG.rds")))
  gene_lists <- get_gene_lists(spe_unweighted, spe_weighted, only_small_ls = T, low = 0.011, high = 0.025)
  
  gene_list1_entrez <- bitr(gene_lists[[1]], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_list2_entrez <- bitr(gene_lists[[2]], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  gene_list <- list(unweighted = gene_list1_entrez$ENTREZID, weighted = gene_list2_entrez$ENTREZID)
  
  compare_ego_bp_small <- compareCluster(geneCluster = gene_list,
                                fun = "enrichGO",
                                OrgDb = org.Hs.eg.db,
                                ont = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.1,
                                readable = TRUE)
  
  compare_ego_cc_small <- compareCluster(geneCluster = gene_list,
                                fun = "enrichGO",
                                OrgDb = org.Hs.eg.db,
                                ont = "CC",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.1,
                                readable = TRUE)
  
  compare_ego_mf_small <- compareCluster(geneCluster = gene_list,
                                fun = "enrichGO",
                                OrgDb = org.Hs.eg.db,
                                ont = "MF",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.1,
                                readable = TRUE)
  
  compare_ego_do_small <- compareCluster(geneCluster = gene_list,
                                fun = "enrichDO",
                                ont = "DO",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.1,
                                readable = TRUE)
  
  gene_lists <- get_gene_lists(spe_unweighted, spe_weighted, only_small_ls = F, low = 0.011, high = 0.025)
  
  gene_list1_entrez <- bitr(gene_lists[[1]], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_list2_entrez <- bitr(gene_lists[[2]], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  gene_list <- list(unweighted = gene_list1_entrez$ENTREZID, weighted = gene_list2_entrez$ENTREZID)
  
  compare_ego_bp_all <- compareCluster(geneCluster = gene_list,
                                fun = "enrichGO",
                                OrgDb = org.Hs.eg.db,
                                ont = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.1,
                                readable = TRUE)
  
  compare_ego_cc_all <- compareCluster(geneCluster = gene_list,
                                fun = "enrichGO",
                                OrgDb = org.Hs.eg.db,
                                ont = "CC",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.1,
                                readable = TRUE)
  
  compare_ego_mf_all <- compareCluster(geneCluster = gene_list,
                                fun = "enrichGO",
                                OrgDb = org.Hs.eg.db,
                                ont = "MF",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.1,
                                readable = TRUE)
  
  compare_ego_do_all <- compareCluster(geneCluster = gene_list,
                                fun = "enrichDO",
                                ont = "DO",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.1,
                                readable = TRUE)
  
  pdf(here("plots", "GO_comparison", paste0("GO_comparison_", d, ".pdf")))
  
  # Visualize the comparison results
  dotplot(compare_ego_bp_all) + ggtitle(paste0(d, " BP Enrichment Comparison (All Genes)"))
  dotplot(compare_ego_bp_small) + ggtitle(paste0(d, " BP Enrichment Comparison (Small Lengthscale)"))
  dotplot(compare_ego_cc_all) + ggtitle(paste0(d, " CC Enrichment Comparison (All Genes)"))
  dotplot(compare_ego_cc_small) + ggtitle(paste0(d, " CC Enrichment Comparison (Small Lengthscale)"))
  dotplot(compare_ego_mf_all) + ggtitle(paste0(d, " MF Enrichment Comparison (All Genes)"))
  dotplot(compare_ego_mf_small) + ggtitle(paste0(d, " MF Enrichment Comparison (Small Lengthscale)"))
  dotplot(compare_ego_do_all) + ggtitle(paste0(d, " DO Enrichment Comparison (All Genes)"))
  dotplot(compare_ego_do_small) + ggtitle(paste0(d, " DO Enrichment Comparison (Small Lengthscale)"))
  
  dev.off()
  
}


# above didn't work so well 
# now let's try just the top 1000 ranked SVGs for weighted vs unweighted

get_gene_lists <- function(spe_unweighted, spe_weighted) {
  
  # get the gene list for the unweighted spe rank <= 1000
  # gene id is in rowData(spe_unweighted)$gene_id
  gene_list_unweighted <- rowData(spe_unweighted)$gene_id[which(rowData(spe_unweighted)$rank <= 1000)]
  
  # get the gene list for the weighted spe using rowData(spe_weighted)$weighted_rank <= 1000
  gene_list_weighted <- rowData(spe_weighted)$gene_id[which(rowData(spe_weighted)$weighted_rank <= 1000)]
  
  # return the gene list for the unweighted and weighted spe
  return(list(gene_list_unweighted, gene_list_weighted))
  
}

dat=c("LobularBreast") 

for (d in dat) {
  print(d)
  spe_unweighted <- readRDS(here("outputs", "results", paste0("spe_human", d, "_nnSVG.rds")))
  spe_weighted <- readRDS(here("outputs", "results", paste0("spe_human", d, "_weighted_nnSVG.rds")))
  gene_lists <- get_gene_lists(spe_unweighted, spe_weighted)
  
  gene_list1_entrez <- bitr(gene_lists[[1]], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_list2_entrez <- bitr(gene_lists[[2]], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  gene_list <- list(unweighted = gene_list1_entrez$ENTREZID, weighted = gene_list2_entrez$ENTREZID)
  
  compare_ego_bp  <- compareCluster(geneCluster = gene_list,
                                         fun = "enrichGO",
                                         OrgDb = org.Hs.eg.db,
                                         ont = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = 0.05,
                                         qvalueCutoff = 0.1,
                                         readable = TRUE)
  
  compare_ego_cc  <- compareCluster(geneCluster = gene_list,
                                         fun = "enrichGO",
                                         OrgDb = org.Hs.eg.db,
                                         ont = "CC",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = 0.05,
                                         qvalueCutoff = 0.1,
                                         readable = TRUE)
  
  compare_ego_mf  <- compareCluster(geneCluster = gene_list,
                                         fun = "enrichGO",
                                         OrgDb = org.Hs.eg.db,
                                         ont = "MF",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = 0.05,
                                         qvalueCutoff = 0.1,
                                         readable = TRUE)
  
  compare_ego_do  <- compareCluster(geneCluster = gene_list,
                                         fun = "enrichDO",
                                         ont = "DO",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = 0.05,
                                         qvalueCutoff = 0.1,
                                         readable = TRUE)
  

  pdf(here("plots", "GO_comparison", paste0("GO_comparison_top1000_", d, ".pdf")))
  
  # Visualize the comparison results
  dotplot(compare_ego_bp ) + ggtitle(paste0(d, " BP Enrichment Comparison "))
  dotplot(compare_ego_cc ) + ggtitle(paste0(d, " CC Enrichment Comparison "))
  dotplot(compare_ego_mf ) + ggtitle(paste0(d, " MF Enrichment Comparison "))
  dotplot(compare_ego_do ) + ggtitle(paste0(d, " DO Enrichment Comparison "))
  
  dev.off()
  
}
