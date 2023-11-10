<h1>  HRDmethyl_public </h1>
BRCA1 and RAD51C promoter methylation exploration in TCGA, ICGC and DepMap datasets

This repo contains streamline processing of TCGA and ICGC. DepMap was processed separately.
** for DepMap preprocessing of rrbs and mRNA, as well as visualisation of methylation status per CpG, scripts are available from scripts/depmap/. details can be found in https://github.com/molonc-lab/rrbs-scripts)
** DepMap figures (Figure 7b, c, d, e) can be found in the plotting.Rmd

# Table of Content
- [Table of Content](#table-of-content)
- [System Requirements](#system-requirements)
  - [R Dependencies](#r-dependencies)
    - [R/4.2.0](#r420)
    - [R/4.0.2 (DepMap analysis)](#r402-depmap-analysis)
    - [R/3.6.2](#r362)
    - [R/3.5.1](#r351)
    - [R/3.4.1](#r341)
  - [Python Dependencies](#python-dependencies)
    - [Python/3.9.13](#python3913)
    - [Python/3.9.6 (DepMap)](#python396-depmap)
- [Required Data](#required-data)
  - [Final Directory Structure](#final-directory-structure)
  - [Raw Data Structure](#raw-data-structure)
    - [IDAT](#idat)
    - [MAF, mRNA](#maf-mrna)
- [running](#running)

# System Requirements
We used CenOS Linux 7 (Core) with
- R/4.2.0, R/3.6.2, R/3.5.1, R/3.4.1 
- Python/3.9.13 Python/3.9.6
- bedtools/2.29.0
- VEP/102
- Bismark/0.19.1

Below are the dependencies for the Python and R languages.

## R Dependencies
### R/4.2.0
- minfi 1.44.0
- data.table 1.14.8
- ggplot2 3.4.2
- tidyverse 2.0.0
- ggnewscale 0.4.9
- edgeR 3.40.2
- org.Hs.eg.db 3.16.0
- scales 1.2.1
- ComplexHeatmap 2.14.0
- rstatix 0.7.2
- ggrepel 0.9.3
- gridExtra 2.3
- readxl 1.4.3
- circlize 0.4.15
- ggbeeswarm 0.7.2
- ggpubr 0.6.0
- purrr 1.0.1

### R/4.0.2 (DepMap analysis)
- ggplot2 3.3.5
- purrr   0.3.4 
- tibble  3.1.6
- dplyr   1.0.10
- tidyr   1.1.4      
- ComplexHeatmap 2.6.2
- methylKit 1.16.1
- edgeR 3.32.1
- RColorBrewer 1.1-2 
- readxl 1.3.1

### R/3.6.2
- ChAMP 2.16.2
- drc 3.0-1

### R/3.5.1
- tidyverse 1.3.1
- BSgenome 1.50.0
- VariantAnnotation 1.28.13
- docopt 0.7.1

### R/3.4.1
- copynumber 1.16.0
- Palimpsest 1.0.0
- R.matlab 3.6.2
- bedr 1.0.4

## Python Dependencies

### Python/3.9.13
- future 0.18.2
- numpy 1.22.4
- pandas 2.0.2
- tqdm 4.65.0
- sklearn 1.1.2
- networkx 2.8.5
- numbpy.linalg
- matplotlib -3.5.3
- seaborn 0.11.2
- scipy 1.9.0
- umap 0.5.3
- netneurotools 0.2.3

### Python/3.9.6 (DepMap)
- matplotlib 3.3.1
- numpy 1.25.2
- pandas 2.0.3


# Required Data
for supporting data from public datasets and publications, please refer to methods and data avaliability section

## Final Directory Structure
```
.
├── data
│   ├── output_error
│   ├── processed
│   │   ├── copynumber
│   │   ├── depmap_brett
│   │   ├── hrdetect_icgc
│   │   ├── hrd_icgc
│   │   │   ├── Array_HRD_score
│   │   │   └── WGS_HRD_score
│   │   ├── methylation
│   │   ├── methylation_clustering
│   │   ├── mRNA
│   │   ├── scarHRD
│   │   ├── Signatures
│   │   └── VCF
│   ├── raw
│   ├── source_data
│   ├── supplementary_data
│   └── supporting
├── figures
├── plots
│   └── Supplementary
└── scripts
```
## Raw Data Structure
raw data availability please refer to the paper
raw data directory structure - example
```
raw
└─── 20220401_TCGA_BLCA_meth
    ├── blca_clinical.tsv
    ├── blca_sample.tsv
    ├── copynumber
    ├── sample.gene_level_copy_number.v36.tsv
    │    └── sample_cnv_sample.tsv
    ├── IDAT
    │    ├── sample_noid_Red.idat
    │    └── sample_noid_Grn.idat
    ├── MAF
    │    ├── sample_wxs.aliquot_ensemble_masked.maf
    │    └── sample_maf_sample.tsv
    └── mRNA
         ├── sample.rna_seq.augmented_star_gene_counts.tsv
         └── sample_mRNA_sample.tsv
```
### IDAT 
IDAT folder contains all the idat readouts from methylation array (HM450 platform)

When checkout the methyltion readouts from GDC portal, blca_clinical.tsv and blca_sample.tsv can be downloaded using the "Sample Sheet" and "Clinical" option.

specifically, the sample sheet is essential for methylation preprocessing.

### MAF, mRNA
apart from masked maf file, a sample sheet for each assay and dataset should be acquired using the method described in the IDAT section

# running
generation of processed data that is required to plot
- step 1. acquire all the raw and supporting data needed
- step 2. for TCGA and ICGC data, run script main.Rmd in the root directory
- step 3. for depmap data, run scripts from scripts/depmap (details please see https://github.com/molonc-lab/rrbs-scripts)
- step 3. run plotting.Rmd in the root directory
