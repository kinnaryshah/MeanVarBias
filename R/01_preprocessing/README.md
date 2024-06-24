# Overview

These are stand alone scripts that allow you to download and preprocess individual data sets use in the `spoon` paper. 

Before running the scripts below, please run the `/R/00_setup.R` script. 

The following datasets are included in this folder: 

- Human dorsolateral preprefrontal cortex (DLPFC)
    - Source: https://research.libd.org/spatialLIBD/
    - Paper: [Maynard et al. 2021](https://www.nature.com/articles/s41593-020-00787-0)
    - Script: `preprocessing_humanDLFPC.R`
- Human breast cancer tissue: 
    - Source: https://www.10xgenomics.com/datasets/human-breast-cancer-visium-fresh-frozen-whole-transcriptome-1-standard
    - Script 1: `preprocessing_humanBreast.sh` (run in cmd to download files from 10x website)
    - Script 2: `preprocessing_humanBreast.R` (run in R to create RDS object)
- Human locus coeruleus (LC): 
    - Source: https://bioconductor.org/packages/WeberDivechaLCdata
    - Paper: [Weber and Divecha et al. 2024](https://elifesciences.org/articles/84628)
    - Script: `preprocessing_humanLC.R`
- Human hippocampus (HPC): 
    - Source: https://github.com/LieberInstitute/spatial_hpc
    - Paper: [Nelson et al. 2024](https://doi.org/10.1101/2024.04.26.590643)
    - Script: `preprocessing_humanHPC.R`
- Human ovarian cancer tissue: 
    - Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211956
        - Tissue sample selected: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6506111
    - Paper: [Denisenko et al. 2024](https://pubmed.ncbi.nlm.nih.gov/38570491/)
    - Script 1: `preprocessing_humanOvarian.sh` (run in cmd to download files from GEO)
    - Script 2: `preprocessing_humanOvarian.R` (run in R to create RDS object)
- Simulated data from a gamma-Poisson distribution (using `splatter`): 
    - Source: https://bioconductor.org/packages/splatter
    - Script: `preprocesing_simulated_Gamma-Poisson.R`

