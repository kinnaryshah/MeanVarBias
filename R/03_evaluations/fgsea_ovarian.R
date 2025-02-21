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

# OVARIAN

spe_nnSVG_Ovarian <- readRDS(here("outputs", "results", "spe_humanOvarian_nnSVG.rds"))
spe_weighted_nnSVG_Ovarian <- readRDS(here("outputs", "results", "spe_humanOvarian_weighted_nnSVG.rds"))

# remove one instance of "HSPA14" for duplicate row names
idx <- which(rowData(spe_weighted_nnSVG_Ovarian)$gene_name == "HSPA14")[2]

spe_nnSVG_Ovarian <- spe_nnSVG_Ovarian[-idx,]
spe_weighted_nnSVG_Ovarian <- spe_weighted_nnSVG_Ovarian[-idx,]

# only keep significant SVGs
spe_nnSVG_Ovarian <- spe_nnSVG_Ovarian[rowData(spe_nnSVG_Ovarian)$padj < .05,]

Ovarian.unweighted.id = mapIds(org.Hs.eg.db, keys=rowData(spe_nnSVG_Ovarian)$gene_name, 
                               keytype="SYMBOL", column="ENTREZID", multiVals = "first")
names(Ovarian.unweighted.id) = rowData(spe_nnSVG_Ovarian)$gene_name
idx <- is.na(Ovarian.unweighted.id)
Ovarian.unweighted.id = Ovarian.unweighted.id[!idx]
pathways.Ovarian.unweighted <- reactomePathways(Ovarian.unweighted.id)
pathways.Ovarian.unweighted <- x.h[names(pathways.Ovarian.unweighted)]

spe_nnSVG_Ovarian <- spe_nnSVG_Ovarian[!idx,]
Ovarian.unweighted.stats = rowData(spe_nnSVG_Ovarian)$padj
names(Ovarian.unweighted.stats) = Ovarian.unweighted.id

set.seed(123)
olig.results.Ovarian.unweighted = fgseaMultilevel(pathways.Ovarian.unweighted, stats=Ovarian.unweighted.stats,
                                                  scoreType="pos", minSize=15, maxSize=500)

olig.results.Ovarian.unweighted$leadingEdge2 = sapply(olig.results.Ovarian.unweighted$leadingEdge, paste, collapse="/")
small_dat <- olig.results.Ovarian.unweighted[,c(1,3)]
print(small_dat[order(small_dat$padj),][c(1:10),])
print(small_dat[order(small_dat$padj),][c(1:3),]$pathway)

# p = 0.09
#1:                                                                        Eukaryotic Translation Elongation

spe_weighted_nnSVG_Ovarian <- spe_weighted_nnSVG_Ovarian[rowData(spe_weighted_nnSVG_Ovarian)$weighted_padj < .05,]

Ovarian.weighted.id = mapIds(org.Hs.eg.db, keys=rowData(spe_weighted_nnSVG_Ovarian)$gene_name, 
                               keytype="SYMBOL", column="ENTREZID", multiVals = "first")
names(Ovarian.weighted.id) = rowData(spe_weighted_nnSVG_Ovarian)$gene_name
idx <- is.na(Ovarian.weighted.id)
Ovarian.weighted.id = Ovarian.weighted.id[!idx]
pathways.Ovarian.weighted <- reactomePathways(Ovarian.weighted.id)
pathways.Ovarian.weighted <- x.h[names(pathways.Ovarian.weighted)]

spe_weighted_nnSVG_Ovarian <- spe_weighted_nnSVG_Ovarian[!idx,]
Ovarian.weighted.stats = rowData(spe_weighted_nnSVG_Ovarian)$weighted_padj
names(Ovarian.weighted.stats) = Ovarian.weighted.id

set.seed(123)
olig.results.Ovarian.weighted = fgseaMultilevel(pathways.Ovarian.weighted, stats=Ovarian.weighted.stats,
                                                  scoreType="pos", minSize=15, maxSize=500)

olig.results.Ovarian.weighted$leadingEdge2 = sapply(olig.results.Ovarian.weighted$leadingEdge, paste, collapse="/")
small_dat <- olig.results.Ovarian.weighted[,c(1,3)]
print(small_dat[order(small_dat$padj),][c(1:10),])
print(small_dat[order(small_dat$padj),][c(1:3),]$pathway)

# p = 5e-6
#  1:                                     SARS-CoV-2 modulates host translation machinery
