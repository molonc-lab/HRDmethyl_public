args = commandArgs(trailingOnly=TRUE)

study_name = "TCGA_CESC"
study_dir = "data/raw/20220401_TCGA_CESC_meth/mRNA"
working_dir = paste0(study_dir, "/../../../../")
plot_output_dir = paste0(working_dir, "plots/mRNA")
output_dir = paste0(working_dir, "data/side_project")

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(edgeR)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))
source(paste0(working_dir,"scripts/mRNA_edgeR_support.R"))

#load methylation file
beta_sub <- fread(paste0(working_dir, "data/processed/methylation/", study_name, "_beta_subset_promo.csv"))
deduplicated_sample_list <- fread(paste0(working_dir, "data/processed/methylation/", study_name, "_deduplicated_sample.csv"), header =T)
beta_sub <- beta_sub %>%
  filter(Sample_Name %in% deduplicated_sample_list$Sample_Name)

#gene annotation
gene_annotation <- fread(paste0(working_dir, "data/supporting/biomaRt_ensg_annotation.csv"), drop =1, header=T)


# prepare count, sample and gene matrix
mRNA_counts <- Sys.glob(paste0(study_dir, "/*gene*ounts*"))%>%
  map_dfr(~as.data.frame(data.table::fread(.x, skip = 6, header = F) %>%
                           mutate(filename = basename(.x)))) %>%
  dplyr::select(gene_ID = V1, gene_name = V2, unstranded_counts = V4, filename)

#sum mRNA count for each file to filter duplicates
sum_counts <- mRNA_counts %>%
  ungroup() %>%
  group_by(filename) %>%
  summarise(sum_counts=sum(unstranded_counts))

gene_csv <- mRNA_counts %>%
  dplyr::select(gene_ID, gene_name) %>%
  mutate(ensembl_gene_id = substr(gene_ID, 1, 15)) %>%
  left_join(gene_annotation, by = "ensembl_gene_id") %>%
  distinct() 

mRNA_counts <- mRNA_counts %>%
  pivot_wider(names_from = filename, values_from = unstranded_counts) %>%
  dplyr::select(-gene_name) %>%
  column_to_rownames("gene_ID")

sample <- as.data.frame(fread(paste0(study_dir, "/", study_name, "_mRNA_sample.tsv")))
#fix column names
colnames(sample) <- make.names(colnames(sample))

#add methylation information for design matrix
sample <- sample %>%
  dplyr::select(filename=File.Name, Sample.ID, Sample.Type)

#filter for higher sum counts among duplicates
filtered_filename <- sum_counts %>%
  left_join(sample, by = c("filename")) %>%
  ungroup() %>%
  group_by(Sample.ID) %>%
  filter(sum_counts == max(sum_counts)) %>%
  ungroup() %>%
  dplyr::select(filename)

sample_list <- sample %>%
  filter(filename %in% unlist(filtered_filename))  #design matrix have BRCA1 and RAD51C_M information combined

#have column matched between mRNA counts and sample information
mRNA_counts <- mRNA_counts[,sample_list$filename]



#form star_matrix
star_matrix <- DGEList(counts = mRNA_counts, group = factor(sample_list$Sample.Type), samples = sample_list)


#QC
##plot raw expression distribution in 50 random sample
ran_sample <- sample(dim(star_matrix)[2], 50)
pdf(paste0(plot_output_dir, "/QC/", study_name, "_raw_distribution50.pdf"))
boxplot(cpm(star_matrix[,ran_sample], log =T), use.cols = T)     
dev.off()

##filter out low count gene
keep <- filterByExpr(star_matrix)
star_matrix <- star_matrix[keep,, keep.lib.sizes=FALSE]
star_matrix <- calcNormFactors(star_matrix)

##re-plot expression distribution
pdf(paste0(plot_output_dir, "/QC/", study_name, "_filter_norm_distribution50.pdf"))
boxplot(cpm(star_matrix[,ran_sample], log =T), use.cols = T)  
dev.off()


# normalisation
design.mat <- model.matrix(~star_matrix$samples$group)
star_matrix <- estimateDisp(star_matrix, design.mat, robust= T)
norm_counts <- cpm(star_matrix, log=TRUE)
fwrite(as.data.table(norm_counts, keep.rownames = ""), paste0(output_dir, "/", study_name, "_mRNA_edgeR.csv"))

##plot BCV
pdf(paste0(plot_output_dir, "/QC/", study_name,"_plotBCV.pdf"))
plotBCV(star_matrix)
dev.off()







