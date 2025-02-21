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

# LobularBreast

spe_nnSVG_LobularBreast <- readRDS(here("outputs", "results", "spe_humanLobularBreast_nnSVG.rds"))
spe_weighted_nnSVG_LobularBreast <- readRDS(here("outputs", "results", "spe_humanLobularBreast_weighted_nnSVG.rds"))

# remove one instance of "HSPA14" for duplicate row names
idx1 <- which(rowData(spe_weighted_nnSVG_LobularBreast)$gene_name == "HSPA14")[2]
idx2 <- which(rowData(spe_weighted_nnSVG_LobularBreast)$gene_name == "CYB561D2")[1]

spe_nnSVG_LobularBreast <- spe_nnSVG_LobularBreast[-c(idx1,idx2),]
spe_weighted_nnSVG_LobularBreast <- spe_weighted_nnSVG_LobularBreast[-c(idx1,idx2),]

# only keep significant SVGs
spe_nnSVG_LobularBreast <- spe_nnSVG_LobularBreast[rowData(spe_nnSVG_LobularBreast)$padj < .05,]

LobularBreast.unweighted.id = mapIds(org.Hs.eg.db, keys=rowData(spe_nnSVG_LobularBreast)$gene_name, 
                                     keytype="SYMBOL", column="ENTREZID", multiVals = "first")
names(LobularBreast.unweighted.id) = rowData(spe_nnSVG_LobularBreast)$gene_name
idx <- is.na(LobularBreast.unweighted.id)
LobularBreast.unweighted.id = LobularBreast.unweighted.id[!idx]
pathways.LobularBreast.unweighted <- reactomePathways(LobularBreast.unweighted.id)
pathways.LobularBreast.unweighted <- x.h[names(pathways.LobularBreast.unweighted)]

spe_nnSVG_LobularBreast <- spe_nnSVG_LobularBreast[!idx,]
LobularBreast.unweighted.stats = rowData(spe_nnSVG_LobularBreast)$padj
names(LobularBreast.unweighted.stats) = LobularBreast.unweighted.id

set.seed(123)
olig.results.LobularBreast.unweighted = fgseaMultilevel(pathways.LobularBreast.unweighted, stats=LobularBreast.unweighted.stats,
                                                        scoreType="pos", minSize=15, maxSize=500)

olig.results.LobularBreast.unweighted$leadingEdge2 = sapply(olig.results.LobularBreast.unweighted$leadingEdge, paste, collapse="/")
small_dat <- olig.results.LobularBreast.unweighted[,c(1,3)]
print(small_dat[order(small_dat$padj),][c(1:10),])
print(small_dat[order(small_dat$padj),][c(1:3),]$pathway)

# even more are sig
#  1:              SRP-dependent cotranslational protein targeting to membrane
#2:                                             Cellular response to hypoxia
#3: Oxygen-dependent proline hydroxylation of Hypoxia-inducible Factor Alpha
#4:                                SCF(Skp2)-mediated degradation of p27/p21
#5:                                    Antigen processing-Cross presentation
#6:                                                     ER-Phagosome pathway
#7:                                              Signaling by ROBO receptors
#8:                                     Cap-dependent Translation Initiation
#9:                                        Eukaryotic Translation Initiation
#10:                                                 Neutrophil degranulation

spe_weighted_nnSVG_LobularBreast <- spe_weighted_nnSVG_LobularBreast[!is.na(rowData(spe_weighted_nnSVG_LobularBreast)$weighted_padj) & rowData(spe_weighted_nnSVG_LobularBreast)$weighted_padj < 0.05, ]

LobularBreast.weighted.id = mapIds(org.Hs.eg.db, keys=rowData(spe_weighted_nnSVG_LobularBreast)$gene_name, 
                                   keytype="SYMBOL", column="ENTREZID", multiVals = "first")
names(LobularBreast.weighted.id) = rowData(spe_weighted_nnSVG_LobularBreast)$gene_name
idx <- is.na(LobularBreast.weighted.id)
LobularBreast.weighted.id = LobularBreast.weighted.id[!idx]
pathways.LobularBreast.weighted <- reactomePathways(LobularBreast.weighted.id)
pathways.LobularBreast.weighted <- x.h[names(pathways.LobularBreast.weighted)]

spe_weighted_nnSVG_LobularBreast <- spe_weighted_nnSVG_LobularBreast[!idx,]
LobularBreast.weighted.stats = rowData(spe_weighted_nnSVG_LobularBreast)$weighted_padj
names(LobularBreast.weighted.stats) = LobularBreast.weighted.id

set.seed(123)
olig.results.LobularBreast.weighted = fgseaMultilevel(pathways.LobularBreast.weighted, stats=LobularBreast.weighted.stats,
                                                      scoreType="pos", minSize=15, maxSize=500)

olig.results.LobularBreast.weighted$leadingEdge2 = sapply(olig.results.LobularBreast.weighted$leadingEdge, paste, collapse="/")
small_dat <- olig.results.LobularBreast.weighted[,c(1,3)]
print(small_dat[order(small_dat$padj),][c(1:10),])
print(small_dat[order(small_dat$padj),][c(1:3),]$pathway)

# even more are significant
#1:                                            Eukaryotic Translation Elongation
#2:                                           Eukaryotic Translation Termination
#3:                                                     Peptide chain elongation
#4:                                                       Viral mRNA Translation
#5:                                                  Selenoamino acid metabolism
#6: Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC)
#7:                  SRP-dependent cotranslational protein targeting to membrane
#8:                                                     Selenocysteine synthesis
#9:                          Response of EIF2AK4 (GCN2) to amino acid deficiency
#10:                                         Cap-dependent Translation Initiation
