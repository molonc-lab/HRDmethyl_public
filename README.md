<h1>  HRDmethyl_public </h1>
BRCA1 and RAD51C promoter methylation exploration in TCGA, ICGC and DepMap datasets


# System requirements
We used CenOS Linux 7 (Core) with
R/4.2.0, R/3.6.2, R/3.5.1, R/3.4.1 
Python/3.9.13
bedtools/2.29.0
VEP/102
Bismark/0.19.1
Below are the dependencies for the Python and R languages.

## R dependencies
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

## Python dependencies
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

# required data
## final directory structure
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

## raw data
raw data availability please refer to the paper
raw data directory structure - example

raw
└─── 20220401_TCGA_BLCA_meth
    ├── copynumber
    ├── IDAT
    ├── MAF
    └── mRNA

# running
generation of processed data that is required to plot
step 1. acquire all the raw and supporting data needed
step 2. run script main.Rmd in the root directory
step 3. run plotting.Rmd in the root directory