args = commandArgs(trailingOnly=TRUE)

study_name = gsub("_meth", "", args[1])
study_dir = paste0(args[2],"/mRNA")
working_dir = paste0(study_dir, "/../../../../")
plot_output_dir = paste0(working_dir, "plots/mRNA")
output_dir = paste0(working_dir, "data/processed/mRNA")

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(edgeR)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))
source(paste0(working_dir,"/scripts/mRNA_edgeR_support.R"))

#load methylation file
beta_sub <- fread(paste0(working_dir, "data/processed/methylation/", study_name, "_beta_subset_promo.csv"))
deduplicated_sample_list <- fread(paste0(working_dir, "data/processed/methylation/", study_name, "_deduplicated_sample.csv"), header =T)
beta_sub <- beta_sub %>%
  filter(Sample_Name %in% deduplicated_sample_list$Sample_Name)

#annotate methylation status for mRNA DEG analysis
beta_sub_meth <- BRCA1_RAD51C_meth_annotat(beta_sub%>%distinct())

#gene annotation
gene_annotation <- fread(paste0(working_dir, "data/supporting/biomaRt_ensg_annotation.csv"), drop =1, header=T)


# prepare count, sample and gene matrix
if(str_detect(study_name, "ICGC")){
  mRNA_counts <- Sys.glob(paste0(study_dir, "/*gene*ounts*"))%>%
    map_dfr(~ as.data.frame(fread(.x))%>%
            dplyr::select(gene_ID = `Ensembl Gene ID`, contains("AOCS")))
  
  gene_csv <- mRNA_counts %>%
    dplyr::select(ensembl_gene_id = gene_ID) %>%
    left_join(gene_annotation, by = "ensembl_gene_id") %>%
    distinct() 
  
  sample <- as.data.frame(fread(paste0(study_dir, "/", study_name, "_mRNA_sample.csv"))) %>%
    dplyr::select(Sample_Name = sampleLabel, Sample.ID = donorLabel, donorLabel, SampleType, sampleDescription) 
  
  sample_list <-as.data.frame(colnames(mRNA_counts)[-1] ) %>%
    mutate(sum = colSums(mRNA_counts %>% dplyr:: select(-gene_ID), na.rm = T)) %>%
    mutate(mRNA_rowname = `colnames(mRNA_counts)[-1]`) %>%
    separate(col = "colnames(mRNA_counts)[-1]", into = c("Sample_Name_mRNA", "Sample.ID"), sep = "_")%>%
    mutate(Sample.ID = substr(Sample.ID, 1, 8)) %>%
    left_join(sample, by =c("Sample_Name_mRNA"="Sample_Name", "Sample.ID")) %>%
    filter(str_detect(SampleType, "Primary")) %>%
    distinct() %>%
    mutate(Sample.ID = gsub("-", "_", substr(Sample.ID, 1, 8)),
           Sample_Name2 =  gsub("-", "_",Sample_Name_mRNA)) %>%
    mutate(Sample_Name_mRNA = paste(Sample.ID, Sample_Name2, sep = "_")) %>%
    dplyr::select(- Sample_Name2) %>%
    left_join(beta_sub_meth %>% 
                mutate(Patient.ID = substr(Sample.ID, 1, 8)) %>%
                filter(Sample_Group == "Primary Tumour") %>%
                dplyr::select(Sample.ID= Patient.ID, BRCA1_M, RAD51C_M) %>%
                distinct())  %>%
    group_by(Sample.ID) %>%
    filter(sum == max(sum)) %>%
    ungroup() %>%
    unite("BRCA1_RAD51C", BRCA1_M:RAD51C_M)
  
  mRNA_counts <- mRNA_counts %>%
    column_to_rownames(var = "gene_ID")
  mRNA_counts <- mRNA_counts[, sample_list$mRNA_rowname]
  
} else {
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
    dplyr::select(filename=File.Name, Sample.ID, Sample.Type) %>%
    mutate(Sample.ID = substr(Sample.ID, 1, 15)) %>%
    left_join(beta_sub_meth %>%
                mutate(Sample.ID = substr(Sample.ID, 1, 15)) %>%
                dplyr::select(Sample.ID, BRCA1_M, RAD51C_M) %>%
                distinct() , by = "Sample.ID") %>%
    distinct() %>%
    mutate(BRCA1_M = ifelse(
      str_detect(Sample.Type, regex('normal', ignore_case = T)), "No", BRCA1_M)) %>%
    mutate(RAD51C_M = ifelse(
      str_detect(Sample.Type, regex('normal', ignore_case = T)), "No", RAD51C_M)) %>%
    filter((!is.na(BRCA1_M)) | !is.na(RAD51C_M))
  
  #filter for higher sum counts among duplicates
  filtered_filename <- sum_counts %>%
    left_join(sample, by = c("filename")) %>%
    ungroup() %>%
    group_by(Sample.ID) %>%
    filter(sum_counts == max(sum_counts)) %>%
    ungroup() %>%
    dplyr::select(filename)
  
  sample_list <- sample %>%
    filter(filename %in% unlist(filtered_filename)) %>%
    unite("BRCA1_RAD51C", BRCA1_M:RAD51C_M) #design matrix have BRCA1 and RAD51C_M information combined
  
  #remove normal
  sample_list <- sample_list %>%
    filter(! str_detect(Sample.Type, regex('normal', ignore_case = T)))
  #have column matched between mRNA counts and sample information
  mRNA_counts <- mRNA_counts[,sample_list$filename]
}


#form star_matrix
star_matrix <- DGEList(counts = mRNA_counts, group = factor(sample_list$BRCA1_RAD51C), samples = sample_list)


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

#DGE
fit <- glmQLFit(star_matrix, design.mat)
if(dim(design.mat)[2] == 3){
  if(sum(design.mat[,2]) >=3){
  et12 <- glmQLFTest(fit, coef=2)
  print("RAD51C regulated")
  print(table(decideTestsDGE(et12, p.value = 0.05)))
  write.csv(topTags(et12, n= 10000, p.value = 0.05),paste0(output_dir, "/", study_name,"_RAD51C_edgeR.csv"), row.names = T)
  print("RAD51C p value")
  print(topTags(et12, n= 10000)["ENSG00000108384",])
  
  #go enrichment
  gene_list <- unlist(mget(substr(rownames(et12),1,15), envir=org.Hs.egENSEMBL2EG,
                           ifnotfound = NA))
  gene_list<- gene_list[substr(rownames(et12),1 ,15)]
  
  go12<- goana(et12, species = "Hs", geneid = gene_list)
  write.csv(go12, paste0(output_dir, "/", study_name,"_RAD51C_go.csv"), row.names = T)
  
  }
  if(sum(design.mat[,3]) >=3){
    et13 <- glmQLFTest(fit, coef=3)
    print("BRCA1 regulated")
    print(table(decideTestsDGE(et13, p.value = 0.05)))
    write.csv(topTags(et13, n= 10000, p.value = 0.05),paste0(output_dir, "/", study_name,"_BRCA1_edgeR.csv"), row.names = T)
    print("BRCA1 p value")
    print(topTags(et13, n= 10000)["ENSG00000012048",])
    
    #go enrichment
    gene_list <- unlist(mget(substr(rownames(et13),1,15), envir=org.Hs.egENSEMBL2EG,
                             ifnotfound = NA))
    gene_list<- gene_list[substr(rownames(et13),1 ,15)]
    
    go13<- goana(et13, species = "Hs", geneid = gene_list)
    write.csv(go13, paste0(output_dir, "/", study_name,"_BRCA1_go.csv"), row.names = T)
  }
}else{
  if(sum(design.mat[,2]) >=3){
    et12 <- glmQLFTest(fit, coef=2)
    print("BRCA1 regulated")
    print(table(decideTestsDGE(et12, p.value = 0.05)))
    write.csv(topTags(et12, n= 10000, p.value = 0.05),paste0(output_dir, "/", study_name,"_BRCA1_edgeR.csv"), row.names = T)
    print("BRCA1 p value")
    print(topTags(et12, n= 10000)["ENSG00000012048",])
    
    #go enrichment
    gene_list <- unlist(mget(substr(rownames(et12),1,15), envir=org.Hs.egENSEMBL2EG,
                             ifnotfound = NA))
    gene_list<- gene_list[substr(rownames(et12),1 ,15)]
    
    go12<- goana(et12, species = "Hs", geneid = gene_list)
    write.csv(go12, paste0(output_dir, "/", study_name,"_BRCA1_go.csv"), row.names = T)
  }
}





