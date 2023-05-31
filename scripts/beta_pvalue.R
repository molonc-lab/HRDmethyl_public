#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

suppressMessages(suppressWarnings(library(ChAMP)))

study_dir = args[1]
idat_dir = paste0(study_dir,"/", "IDAT")
myLoad <- champ.load(directory = idat_dir, method = "ChaAMP", arraytype = "450K")


