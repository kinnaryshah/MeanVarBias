library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(here)
library(SpatialExperiment)
library(ggplot2)
library(enrichplot)
library(patchwork)

# create a function to get a list of genes from the argument spe object
# we want to get the list of genes that are SVGs for the weighted and unweighted versions of the spe object
# if only_small_ls is TRUE, we want the small lengthscale genes only
# if only_small_ls is FALSE, we want all genes
# return a list of genes for the unweighted and weighted versions of the spe object that can be input in GO analysis

get_gene_lists <- function(spe_unweighted, spe_weighted, p_val_thres = F, n_SVGs = 1000) {
  if (p_val_thres) {
    # get the gene list for the unweighted spe using rowData(spe_unweighted)$rank < 1000
    idx <- which(rowData(spe_unweighted)$rank <= n_SVGs)
    gene_list_unweighted <- rowData(spe_unweighted)$LR_stat[idx]
    names(gene_list_unweighted) <- rowData(spe_unweighted)$gene_id[idx]
    
    # get the gene list for the weighted spe using rowData(spe_weighted)$weighted_rank < 1000
    idx <- which(rowData(spe_weighted)$weighted_rank <= n_SVGs)
    gene_list_weighted <- rowData(spe_weighted)$weighted_LR_stat[idx]
    names(gene_list_weighted) <- rowData(spe_weighted)$gene_id[idx]
    
  }
  
  else {
    # get the gene list for the unweighted spe using rowData(spe_unweighted)$padj < 0.05
    idx <- which(rowData(spe_unweighted)$padj < 0.05)
    gene_list_unweighted <- rowData(spe_unweighted)$LR_stat[idx]
    names(gene_list_unweighted) <- rowData(spe_unweighted)$gene_id[idx]
    
    # get the gene list for the weighted spe using rowData(spe_weighted)$weighted_padj < 0.05
    idx <- which(rowData(spe_weighted)$weighted_padj < 0.05)
    gene_list_weighted <- rowData(spe_weighted)$weighted_LR_stat[idx]
    names(gene_list_weighted) <- rowData(spe_weighted)$gene_id[idx]

  }
  
  # return the gene list for the unweighted and weighted spe
  return(list(gene_list_unweighted, gene_list_weighted))
  
}

dat=c("LobularBreast", "Breast", "Ovarian", "SubtypeBreast")
plot_list <- list()

for (d in dat) {
  print(d)
  spe_unweighted <- readRDS(here("outputs", "results", paste0("spe_human", d, "_nnSVG.rds")))
  spe_weighted <- readRDS(here("outputs", "results", paste0("spe_human", d, "_weighted_nnSVG.rds")))
  
  gene_lists <- get_gene_lists(spe_unweighted, spe_weighted, p_val_thres = T, n_SVGs = 1000)
  
  gene_list1_entrez <- bitr(names(gene_lists[[1]]), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene_list2_entrez <- bitr(names(gene_lists[[2]]), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  ensembl_to_entrez1 <- setNames(gene_list1_entrez$ENTREZID, gene_list1_entrez$ENSEMBL)
  
  mapped_names1 <- ensembl_to_entrez1[names(gene_lists[[1]])]
  valid_idx1 <- !is.na(mapped_names1)
  gene_lists[[1]] <- gene_lists[[1]][valid_idx1] 
  names(gene_lists[[1]]) <- mapped_names1[valid_idx1]
  
  ensembl_to_entrez2 <- setNames(gene_list2_entrez$ENTREZID, gene_list2_entrez$ENSEMBL)
  
  mapped_names2 <- ensembl_to_entrez2[names(gene_lists[[2]])]
  valid_idx2 <- !is.na(mapped_names2)
  gene_lists[[2]] <- gene_lists[[2]][valid_idx2]
  names(gene_lists[[2]]) <- mapped_names2[valid_idx2]
  
  # remove duplicated genes
  gene_lists[[1]] <- gene_lists[[1]][!duplicated(names(gene_lists[[1]]) )]
  gene_lists[[2]] <- gene_lists[[2]][!duplicated(names(gene_lists[[2]]) )]
  
  gene_list <- list(unweighted = gene_lists[[1]], weighted = gene_lists[[2]])
  
  set.seed(123)
  gseDGN_small_unweighted <- gseDGN(sort(gene_list[[1]], decreasing = T),
                                  minGSSize     = 15,
                                  pvalueCutoff  = 0.25,
                                  pAdjustMethod = "BH",
                                  scoreType="pos",
                                  verbose       = T)
  gseDGN_small_unweighted <- setReadable(gseDGN_small_unweighted, 'org.Hs.eg.db')
  head(gseDGN_small_unweighted,10)
  
  set.seed(123)
  gseDGN_small_weighted <- gseDGN(sort(gene_list[[2]], decreasing = T),
                                  minGSSize     = 15,
                                  pvalueCutoff  = 0.25,
                                  pAdjustMethod = "BH",
                                  scoreType="pos",
                                  verbose       = T)
  gseDGN_small_weighted <- setReadable(gseDGN_small_weighted, 'org.Hs.eg.db')
  head(gseDGN_small_weighted,10)

  title_d <- d
  
  # Visualize the comparison results
  if (d=="Breast"){
    title_d <- "Ductal Breast"
  } 
  
  if (d=="SubtypeBreast"){
    title_d <- "ER+ Breast"
  } 
  
  if (d=="LobularBreast"){
    title_d <- "Lobular Breast"
  } 
  
  p1 <- dotplot(gseDGN_small_unweighted) + ggtitle(paste0(title_d, " Unweighted nnSVG GSEA"))
  p2 <- dotplot(gseDGN_small_weighted) + ggtitle(paste0(title_d, " Weighted nnSVG GSEA"))
  
  plot_list[[d]] <- list(p1, p2)
  
}

p <- wrap_plots(list(
  plot_list[["LobularBreast"]][[1]],
  plot_list[["LobularBreast"]][[2]],
  plot_list[["Breast"]][[1]],
  plot_list[["Breast"]][[2]],
  plot_list[["Ovarian"]][[1]],
  plot_list[["Ovarian"]][[2]],
  plot_list[["SubtypeBreast"]][[1]],
  plot_list[["SubtypeBreast"]][[2]]
), nrow = 4)

pdf(here("plots", "GO_comparison", "GO_comparison_full.pdf"), height = 16, width = 14)
print(p)
dev.off()
