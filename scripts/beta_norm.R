#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

suppressMessages(suppressWarnings(library(minfi)))
suppressMessages(suppressWarnings(library(data.table)))

study_name = args[1]
study_dir = args[2]

cancer_type = unlist(strsplit(study_name,"_"))[2]
idat_dir = paste0(study_dir,"/", "IDAT")

targets <- read.metharray.sheet(idat_dir)
RGset <- read.metharray.exp(targets = targets)
sampleNames(RGset) <- gsub("_noid", "", sampleNames(RGset))
grSet <- preprocessFunnorm(RGset)
norm_beta <- getBeta(grSet)
fwrite(as.data.table(norm_beta, keep.rownames = ""), paste0("data/processed/methylation/", study_name, "_norm_beta.csv"), row.names = F, sep = "\t")


