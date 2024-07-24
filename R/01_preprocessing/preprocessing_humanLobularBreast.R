###################################
# dataset: Visium human breast cancer 
# https://www.10xgenomics.com/datasets/human-breast-cancer-whole-transcriptome-analysis-1-standard-1-2-0
###################################
library(scran)
library(SpatialExperiment)
library(here)
library(nnSVG)
# BiocManager::install("DropletUtils") if not installed


# ---------
# load data
# ---------

spe <-  read10xVisium(here("outputs/raw/humanLobularBreast/outs"),
                      type = "sparse",   # use sparse (not HDF5) format
                      data = "raw",     
                      images = "lowres", # specify which image(s) to include
                      load = TRUE)      # specify whether or not to load image(s)

# -------------
# preprocessing
# -------------

#row data formatting
rowData(spe)$gene_id <- rownames(spe)
rowData(spe)$gene_name <- rowData(spe)$symbol
rowData(spe)$symbol <- NULL

#remove 667 spots not in tissue
spe <- spe[, colData(spe)$in_tissue == 1]

## remove genes without enough data
spe <- filter_genes(spe, 2, 0.2)
dim(spe)
# > dim(spe)
# [1] 12624  4325

# create logcounts matrix
spe <- logNormCounts(spe)


# -----------
# save object
# -----------

fn <- here("outputs", "preprocessed", "spe_humanLobularBreast_preprocessed.rds")
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
#   [1] splatter_1.29.0             WeberDivechaLCdata_1.7.0   
# [3] ExperimentHub_2.13.0        AnnotationHub_3.13.0       
# [5] BiocFileCache_2.13.0        dbplyr_2.5.0               
# [7] ggspavis_1.11.0             ggplot2_3.5.1              
# [9] here_1.0.1                  nnSVG_1.9.0                
# [11] spatialLIBD_1.17.2          scran_1.33.0               
# [13] scuttle_1.15.0              SpatialExperiment_1.15.1   
# [15] SingleCellExperiment_1.27.2 SummarizedExperiment_1.35.0
# [17] Biobase_2.65.0              GenomicRanges_1.57.1       
# [19] GenomeInfoDb_1.41.1         IRanges_2.39.0             
# [21] S4Vectors_0.43.0            BiocGenerics_0.51.0        
# [23] MatrixGenerics_1.17.0       matrixStats_1.3.0          
# 
# loaded via a namespace (and not attached):
#   [1] later_1.3.2               BiocIO_1.15.0            
# [3] bitops_1.0-7              filelock_1.0.3           
# [5] fields_15.2               R.oo_1.26.0              
# [7] tibble_3.2.1              XML_3.99-0.16.1          
# [9] lifecycle_1.0.4           rprojroot_2.0.4          
# [11] edgeR_4.3.4               doParallel_1.0.17        
# [13] lattice_0.22-6            backports_1.5.0          
# [15] magrittr_2.0.3            limma_3.61.2             
# [17] plotly_4.10.4             sass_0.4.9               
# [19] jquerylib_0.1.4           yaml_2.3.8               
# [21] metapod_1.13.0            httpuv_1.6.15            
# [23] ggside_0.3.1              spam_2.10-0              
# [25] sessioninfo_1.2.2         pbapply_1.7-2            
# [27] cowplot_1.1.3             DBI_1.2.3                
# [29] RColorBrewer_1.1-3        golem_0.4.1              
# [31] maps_3.4.2                abind_1.4-5              
# [33] zlibbioc_1.51.1           R.utils_2.12.3           
# [35] purrr_1.0.2               RCurl_1.98-1.14          
# [37] rappdirs_0.3.3            GenomeInfoDbData_1.2.12  
# [39] ggrepel_0.9.5             irlba_2.3.5.1            
# [41] dqrng_0.4.1               DelayedMatrixStats_1.27.1
# [43] DropletUtils_1.25.0       codetools_0.2-20         
# [45] DelayedArray_0.31.3       DT_0.33                  
# [47] tidyselect_1.2.1          UCSC.utils_1.1.0         
# [49] ScaledMatrix_1.13.0       viridis_0.6.5            
# [51] shinyWidgets_0.8.6        GenomicAlignments_1.41.0 
# [53] jsonlite_1.8.8            BiocNeighbors_1.23.0     
# [55] scater_1.33.1             iterators_1.0.14         
# [57] foreach_1.5.2             tools_4.4.1              
# [59] Rcpp_1.0.12               glue_1.7.0               
# [61] gridExtra_2.3             SparseArray_1.5.9        
# [63] HDF5Array_1.33.3          dplyr_1.1.4              
# [65] withr_3.0.0               BiocManager_1.30.23      
# [67] fastmap_1.2.0             rhdf5filters_1.17.0      
# [69] bluster_1.15.0            fansi_1.0.6              
# [71] digest_0.6.35             rsvd_1.0.5               
# [73] R6_2.5.1                  mime_0.12                
# [75] colorspace_2.1-0          RSQLite_2.3.7            
# [77] R.methodsS3_1.8.2         config_0.3.2             
# [79] utf8_1.2.4                tidyr_1.3.1              
# [81] generics_0.1.3            data.table_1.15.4        
# [83] rtracklayer_1.65.0        httr_1.4.7               
# [85] htmlwidgets_1.6.4         S4Arrays_1.5.1           
# [87] pkgconfig_2.0.3           gtable_0.3.5             
# [89] rdist_0.0.5               blob_1.2.4               
# [91] BRISC_1.0.5               XVector_0.45.0           
# [93] htmltools_0.5.8.1         dotCall64_1.1-1          
# [95] scales_1.3.0              png_0.1-8                
# [97] attempt_0.3.1             rstudioapi_0.16.0        
# [99] rjson_0.2.21              checkmate_2.3.1          
# [101] curl_5.2.1                rhdf5_2.49.0             
# [103] cachem_1.1.0              stringr_1.5.1            
# [105] BiocVersion_3.20.0        parallel_4.4.1           
# [107] vipor_0.4.7               AnnotationDbi_1.67.0     
# [109] restfulr_0.0.15           pillar_1.9.0             
# [111] grid_4.4.1                vctrs_0.6.5              
# [113] RANN_2.6.1                promises_1.3.0           
# [115] BiocSingular_1.21.1       beachmat_2.21.3          
# [117] xtable_1.8-4              cluster_2.1.6            
# [119] beeswarm_0.4.0            paletteer_1.6.0          
# [121] magick_2.8.3              cli_3.6.3                
# [123] locfit_1.5-9.9            compiler_4.4.1           
# [125] Rsamtools_2.21.0          rlang_1.1.4              
# [127] crayon_1.5.3              rematch2_2.1.2           
# [129] ggbeeswarm_0.7.2          stringi_1.8.4            
# [131] viridisLite_0.4.2         BiocParallel_1.39.0      
# [133] munsell_0.5.1             Biostrings_2.73.1        
# [135] lazyeval_0.2.2            Matrix_1.7-0             
# [137] benchmarkme_1.0.8         sparseMatrixStats_1.17.2 
# [139] bit64_4.0.5               Rhdf5lib_1.27.0          
# [141] KEGGREST_1.45.1           statmod_1.5.0            
# [143] shiny_1.8.1.1             igraph_2.0.3             
# [145] memoise_2.0.1             bslib_0.7.0              
# [147] benchmarkmeData_1.0.4     bit_4.0.5    
