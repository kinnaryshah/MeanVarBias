###################################
# dataset: Visium human locus coeruleus
###################################
library(SpatialExperiment)
library(ggspavis)
library(WeberDivechaLCdata)
library(here)

# ---------
# load data
# ---------
spe <- WeberDivechaLCdata_Visium()
#subset one sample
spe <- spe[, spe$sample_id == "Br2701_LC_round2"]
dim(spe)
# > dim(spe)
# [1] 23728  2809


# -------------
# preprocessing
# https://github.com/lmweber/locus-c/blob/main/code/analyses/Visium/06_LC_gene_filtering.R
# -------------

# filter out genes with extremely low expression
# using simple threshold on total UMI counts summed across all spots
n_umis <- 80
ix_low_genes <- rowSums(counts(spe)) < n_umis
table(ix_low_genes)

spe <- spe[!ix_low_genes, ]
dim(spe)
# > dim(spe)
# [1] 1331 2809

# -----------
# save object
# -----------

fn <- here("outputs", "preprocessed", "spe_humanLC_preprocessed.rds")
saveRDS(spe, file = fn)

# -----------
# session information
# -----------

sessionInfo()

# > sessionInfo()
# R version 4.4.1 Patched (2024-06-20 r86815)
# Platform: x86_64-conda-linux-gnu
# Running under: Rocky Linux 9.2 (Blue Onyx)
# 
# Matrix products: default
# BLAS:   /jhpce/shared/community/core/conda_R/4.4.x/R/lib64/R/lib/libRblas.so 
# LAPACK: /jhpce/shared/community/core/conda_R/4.4.x/R/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: US/Eastern
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices datasets  utils     methods  
# [8] base     
# 
# other attached packages:
#   [1] WeberDivechaLCdata_1.7.0    ExperimentHub_2.13.0       
# [3] AnnotationHub_3.13.0        BiocFileCache_2.13.0       
# [5] dbplyr_2.5.0                ggspavis_1.11.0            
# [7] ggplot2_3.5.1               here_1.0.1                 
# [9] nnSVG_1.9.0                 spatialLIBD_1.17.2         
# [11] scran_1.33.0                scuttle_1.15.0             
# [13] SpatialExperiment_1.15.1    SingleCellExperiment_1.27.2
# [15] SummarizedExperiment_1.35.0 Biobase_2.65.0             
# [17] GenomicRanges_1.57.1        GenomeInfoDb_1.41.1        
# [19] IRanges_2.39.0              S4Vectors_0.43.0           
# [21] BiocGenerics_0.51.0         MatrixGenerics_1.17.0      
# [23] matrixStats_1.3.0          
# 
# loaded via a namespace (and not attached):
#   [1] later_1.3.2               BiocIO_1.15.0            
# [3] bitops_1.0-7              filelock_1.0.3           
# [5] fields_15.2               tibble_3.2.1             
# [7] XML_3.99-0.16.1           lifecycle_1.0.4          
# [9] rprojroot_2.0.4           edgeR_4.3.4              
# [11] doParallel_1.0.17         lattice_0.22-6           
# [13] magrittr_2.0.3            limma_3.61.2             
# [15] plotly_4.10.4             sass_0.4.9               
# [17] jquerylib_0.1.4           yaml_2.3.8               
# [19] metapod_1.13.0            httpuv_1.6.15            
# [21] ggside_0.3.1              spam_2.10-0              
# [23] sessioninfo_1.2.2         pbapply_1.7-2            
# [25] cowplot_1.1.3             DBI_1.2.3                
# [27] RColorBrewer_1.1-3        golem_0.4.1              
# [29] maps_3.4.2                abind_1.4-5              
# [31] zlibbioc_1.51.1           purrr_1.0.2              
# [33] RCurl_1.98-1.14           rappdirs_0.3.3           
# [35] GenomeInfoDbData_1.2.12   ggrepel_0.9.5            
# [37] irlba_2.3.5.1             dqrng_0.4.1              
# [39] DelayedMatrixStats_1.27.1 codetools_0.2-20         
# [41] DelayedArray_0.31.3       DT_0.33                  
# [43] tidyselect_1.2.1          UCSC.utils_1.1.0         
# [45] ScaledMatrix_1.13.0       viridis_0.6.5            
# [47] shinyWidgets_0.8.6        GenomicAlignments_1.41.0 
# [49] jsonlite_1.8.8            BiocNeighbors_1.23.0     
# [51] scater_1.33.1             iterators_1.0.14         
# [53] foreach_1.5.2             tools_4.4.1              
# [55] Rcpp_1.0.12               glue_1.7.0               
# [57] gridExtra_2.3             SparseArray_1.5.9        
# [59] dplyr_1.1.4               withr_3.0.0              
# [61] BiocManager_1.30.23       fastmap_1.2.0            
# [63] bluster_1.15.0            fansi_1.0.6              
# [65] digest_0.6.35             rsvd_1.0.5               
# [67] R6_2.5.1                  mime_0.12                
# [69] colorspace_2.1-0          RSQLite_2.3.7            
# [71] config_0.3.2              utf8_1.2.4               
# [73] tidyr_1.3.1               generics_0.1.3           
# [75] data.table_1.15.4         rtracklayer_1.65.0       
# [77] httr_1.4.7                htmlwidgets_1.6.4        
# [79] S4Arrays_1.5.1            pkgconfig_2.0.3          
# [81] gtable_0.3.5              rdist_0.0.5              
# [83] blob_1.2.4                BRISC_1.0.5              
# [85] XVector_0.45.0            htmltools_0.5.8.1        
# [87] dotCall64_1.1-1           scales_1.3.0             
# [89] png_0.1-8                 attempt_0.3.1            
# [91] rjson_0.2.21              curl_5.2.1               
# [93] cachem_1.1.0              stringr_1.5.1            
# [95] BiocVersion_3.20.0        parallel_4.4.1           
# [97] vipor_0.4.7               AnnotationDbi_1.67.0     
# [99] restfulr_0.0.15           pillar_1.9.0             
# [101] grid_4.4.1                vctrs_0.6.5              
# [103] RANN_2.6.1                promises_1.3.0           
# [105] BiocSingular_1.21.1       beachmat_2.21.3          
# [107] xtable_1.8-4              cluster_2.1.6            
# [109] beeswarm_0.4.0            paletteer_1.6.0          
# [111] magick_2.8.3              cli_3.6.3                
# [113] locfit_1.5-9.9            compiler_4.4.1           
# [115] Rsamtools_2.21.0          rlang_1.1.4              
# [117] crayon_1.5.3              rematch2_2.1.2           
# [119] ggbeeswarm_0.7.2          stringi_1.8.4            
# [121] viridisLite_0.4.2         BiocParallel_1.39.0      
# [123] munsell_0.5.1             Biostrings_2.73.1        
# [125] lazyeval_0.2.2            Matrix_1.7-0             
# [127] benchmarkme_1.0.8         sparseMatrixStats_1.17.2 
# [129] bit64_4.0.5               KEGGREST_1.45.1          
# [131] statmod_1.5.0             shiny_1.8.1.1            
# [133] igraph_2.0.3              memoise_2.0.1            
# [135] bslib_0.7.0               benchmarkmeData_1.0.4    
# [137] bit_4.0.5                

