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

# SubtypeBreast

spe_nnSVG_SubtypeBreast <- readRDS(here("outputs", "results", "spe_humanSubtypeBreast_nnSVG.rds"))
spe_weighted_nnSVG_SubtypeBreast <- readRDS(here("outputs", "results", "spe_humanSubtypeBreast_weighted_nnSVG.rds"))

# remove one instance of "HSPA14" for duplicate row names
idx1 <- which(rowData(spe_weighted_nnSVG_SubtypeBreast)$gene_name == "HSPA14")[2]
idx2 <- which(rowData(spe_weighted_nnSVG_SubtypeBreast)$gene_name == "CYB561D2")[1]

spe_nnSVG_SubtypeBreast <- spe_nnSVG_SubtypeBreast[-c(idx1,idx2),]
spe_weighted_nnSVG_SubtypeBreast <- spe_weighted_nnSVG_SubtypeBreast[-c(idx1,idx2),]

# only keep significant SVGs
spe_nnSVG_SubtypeBreast <- spe_nnSVG_SubtypeBreast[rowData(spe_nnSVG_SubtypeBreast)$padj < .05,]

SubtypeBreast.unweighted.id = mapIds(org.Hs.eg.db, keys=rowData(spe_nnSVG_SubtypeBreast)$gene_name, 
                               keytype="SYMBOL", column="ENTREZID", multiVals = "first")
names(SubtypeBreast.unweighted.id) = rowData(spe_nnSVG_SubtypeBreast)$gene_name
idx <- is.na(SubtypeBreast.unweighted.id)
SubtypeBreast.unweighted.id = SubtypeBreast.unweighted.id[!idx]
pathways.SubtypeBreast.unweighted <- reactomePathways(SubtypeBreast.unweighted.id)
pathways.SubtypeBreast.unweighted <- x.h[names(pathways.SubtypeBreast.unweighted)]

spe_nnSVG_SubtypeBreast <- spe_nnSVG_SubtypeBreast[!idx,]
SubtypeBreast.unweighted.stats = rowData(spe_nnSVG_SubtypeBreast)$padj
names(SubtypeBreast.unweighted.stats) = SubtypeBreast.unweighted.id

set.seed(123)
olig.results.SubtypeBreast.unweighted = fgseaMultilevel(pathways.SubtypeBreast.unweighted, stats=SubtypeBreast.unweighted.stats,
                                                  scoreType="pos", minSize=15, maxSize=500)

olig.results.SubtypeBreast.unweighted$leadingEdge2 = sapply(olig.results.SubtypeBreast.unweighted$leadingEdge, paste, collapse="/")
small_dat <- olig.results.SubtypeBreast.unweighted[,c(1,3)]
print(small_dat[order(small_dat$padj),][c(1:10),])
print(small_dat[order(small_dat$padj),][c(1:3),]$pathway)

# even more are sig
# 1:                                                       Influenza Infection
#2:                         Influenza Viral RNA Transcription and Replication
#3:                       Response of EIF2AK4 (GCN2) to amino acid deficiency
#4:                                               Selenoamino acid metabolism
#5:                               Regulation of expression of SLITs and ROBOs
#6: Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC)
#7:                                             Nonsense-Mediated Decay (NMD)
#8:                                                Mitotic Spindle Checkpoint
#9:                                               Arachidonic acid metabolism
#10:                          Cyclin A:Cdk2-associated events at S phase entry

spe_weighted_nnSVG_SubtypeBreast <- spe_weighted_nnSVG_SubtypeBreast[!is.na(rowData(spe_weighted_nnSVG_SubtypeBreast)$weighted_padj) & rowData(spe_weighted_nnSVG_SubtypeBreast)$weighted_padj < 0.05, ]

SubtypeBreast.weighted.id = mapIds(org.Hs.eg.db, keys=rowData(spe_weighted_nnSVG_SubtypeBreast)$gene_name, 
                             keytype="SYMBOL", column="ENTREZID", multiVals = "first")
names(SubtypeBreast.weighted.id) = rowData(spe_weighted_nnSVG_SubtypeBreast)$gene_name
idx <- is.na(SubtypeBreast.weighted.id)
SubtypeBreast.weighted.id = SubtypeBreast.weighted.id[!idx]
pathways.SubtypeBreast.weighted <- reactomePathways(SubtypeBreast.weighted.id)
pathways.SubtypeBreast.weighted <- x.h[names(pathways.SubtypeBreast.weighted)]

spe_weighted_nnSVG_SubtypeBreast <- spe_weighted_nnSVG_SubtypeBreast[!idx,]
SubtypeBreast.weighted.stats = rowData(spe_weighted_nnSVG_SubtypeBreast)$weighted_padj
names(SubtypeBreast.weighted.stats) = SubtypeBreast.weighted.id

set.seed(123)
olig.results.SubtypeBreast.weighted = fgseaMultilevel(pathways.SubtypeBreast.weighted, stats=SubtypeBreast.weighted.stats,
                                                scoreType="pos", minSize=15, maxSize=500)

olig.results.SubtypeBreast.weighted$leadingEdge2 = sapply(olig.results.SubtypeBreast.weighted$leadingEdge, paste, collapse="/")
small_dat <- olig.results.SubtypeBreast.weighted[,c(1,3)]
print(small_dat[order(small_dat$padj),][c(1:10),])
print(small_dat[order(small_dat$padj),][c(1:3),]$pathway)

# none are significant