library(ggplot2)
library(SpatialExperiment)
library(patchwork)
library(dplyr)
library(here)
library(fgsea)
library(reactome.db)
library(org.Hs.eg.db)

# prep reactome
xx <- as.list(reactomePATHID2NAME)
reactome.h <- xx[grep("^Homo",xx)]
x <- as.list(reactomePATHID2EXTID)
reactome.h = reactome.h[intersect(names(reactome.h), names(x))]
x.h <- x[names(reactome.h)]
identical(names(x.h), names(reactome.h))
reactome.h = gsub("Homo sapiens: ","",reactome.h)
names(x.h) = reactome.h

# Breast

spe_nnSVG_Breast <- readRDS(here("outputs", "results", "spe_humanBreast_nnSVG.rds"))
spe_weighted_nnSVG_Breast <- readRDS(here("outputs", "results", "spe_humanBreast_weighted_nnSVG.rds"))

# remove one instance of "HSPA14" for duplicate row names
idx1 <- which(rowData(spe_weighted_nnSVG_Breast)$gene_name == "HSPA14")[2]
idx2 <- which(rowData(spe_weighted_nnSVG_Breast)$gene_name == "CYB561D2")[1]

spe_nnSVG_LobularBreast <- spe_nnSVG_Breast[-c(idx1,idx2),]
spe_weighted_nnSVG_LobularBreast <- spe_weighted_nnSVG_Breast[-c(idx1,idx2),]

# only keep significant SVGs
spe_nnSVG_Breast <- spe_nnSVG_Breast[rowData(spe_nnSVG_Breast)$padj < .05,]

Breast.unweighted.id = mapIds(org.Hs.eg.db, keys=rowData(spe_nnSVG_Breast)$gene_name, 
                               keytype="SYMBOL", column="ENTREZID", multiVals = "first")
names(Breast.unweighted.id) = rowData(spe_nnSVG_Breast)$gene_name
idx <- is.na(Breast.unweighted.id)
Breast.unweighted.id = Breast.unweighted.id[!idx]
pathways.Breast.unweighted <- reactomePathways(Breast.unweighted.id)
pathways.Breast.unweighted <- x.h[names(pathways.Breast.unweighted)]

spe_nnSVG_Breast <- spe_nnSVG_Breast[!idx,]
Breast.unweighted.stats = rowData(spe_nnSVG_Breast)$LR_stat
names(Breast.unweighted.stats) = Breast.unweighted.id

set.seed(123)
olig.results.Breast.unweighted = fgseaMultilevel(pathways.Breast.unweighted, stats=Breast.unweighted.stats,
                                                  scoreType="pos", minSize=15, maxSize=500)

olig.results.Breast.unweighted$leadingEdge2 = sapply(olig.results.Breast.unweighted$leadingEdge, paste, collapse="/")
small_dat <- olig.results.Breast.unweighted[,c(1,3)]
print(small_dat[order(small_dat$padj),][c(1:10),])
print(small_dat[order(small_dat$padj),][c(1:3),]$pathway)

spe_weighted_nnSVG_Breast <- spe_weighted_nnSVG_Breast[!is.na(rowData(spe_weighted_nnSVG_Breast)$weighted_padj) & rowData(spe_weighted_nnSVG_Breast)$weighted_padj < .05,]

Breast.weighted.id = mapIds(org.Hs.eg.db, keys=rowData(spe_weighted_nnSVG_Breast)$gene_name, 
                             keytype="SYMBOL", column="ENTREZID", multiVals = "first")
names(Breast.weighted.id) = rowData(spe_weighted_nnSVG_Breast)$gene_name
idx <- is.na(Breast.weighted.id)
Breast.weighted.id = Breast.weighted.id[!idx]
pathways.Breast.weighted <- reactomePathways(Breast.weighted.id)
pathways.Breast.weighted <- x.h[names(pathways.Breast.weighted)]

spe_weighted_nnSVG_Breast <- spe_weighted_nnSVG_Breast[!idx,]
Breast.weighted.stats = rowData(spe_weighted_nnSVG_Breast)$weighted_LR_stat
names(Breast.weighted.stats) = Breast.weighted.id

set.seed(123)
olig.results.Breast.weighted = fgseaMultilevel(pathways.Breast.weighted, stats=Breast.weighted.stats,
                                                scoreType="pos", minSize=15, maxSize=500)

olig.results.Breast.weighted$leadingEdge2 = sapply(olig.results.Breast.weighted$leadingEdge, paste, collapse="/")
small_dat <- olig.results.Breast.weighted[,c(1,3)]
print(small_dat[order(small_dat$padj),][c(1:10),])
print(small_dat[order(small_dat$padj),][c(1:3),]$pathway)
