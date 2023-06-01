#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

suppressMessages(suppressWarnings(library(minfi)))
suppressMessages(suppressWarnings(library(data.table)))

study_name = args[1]
study_dir = args[2]

cancer_type = unlist(strsplit(study_name,"_"))[2]
idat_dir = paste0(study_dir,"/", "IDAT")

# load IDAT
targets <- read.metharray.sheet(idat_dir)
# read intensity
RGset <- read.metharray.exp(targets = targets)
sampleNames(RGset) <- gsub("_noid", "", sampleNames(RGset))
# normalise the intensity using functional normalisation
grSet <- preprocessFunnorm(RGset)
#remove SNP related CpG
grSet <- dropLociWithSnps(grSet, snps=c("SBE","CpG"), maf=0)
# get beta
norm_beta <- getBeta(grSet)
# save file
fwrite(as.data.table(norm_beta, keep.rownames = ""), paste0("data/processed/methylation/", study_name, "_norm_beta.csv"), row.names = F, sep = "\t")


