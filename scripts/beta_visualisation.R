args = commandArgs(trailingOnly=TRUE)

study_name = gsub("_meth", "", args[1])

directory = args[2]

genecodefile = paste0(directory, "/data/supporting/HM450.hg38.manifest.gencode.v36.tsv.gz")
gene_promoter_info = paste0(directory, "/data/supporting/EPDnew_human_ver006_hg38.tsv")
HR_gene_list = c("BRCA1","BRCA2","RAD51C","RAD51D","BRIP1","PALB2","CHEK2")

# load package
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(require(ggnewscale)))
#load supporting functions
source(paste0(directory,"/scripts/beta_visualisation_support.R"))

#load beta
beta_matrix <- Sys.glob(paste0(directory, "/data/processed/methylation/", study_name, "_meth_norm_beta.csv")) %>%
  map_dfr(~ as.data.frame(fread(.x))) %>%
  mutate(CpG = V1)%>%
  column_to_rownames("V1")

# load sample
sample_sheet <- as.data.frame(fread(Sys.glob(paste0(directory, "/data/raw/2022*_", study_name, "_meth/IDAT/*_pd.csv")), header = T))
# load pvalue
pvalue <- as.data.frame(fread(Sys.glob(paste0(directory, "/data/processed/methylation/", study_name, "_beta_pvalue.csv")), header = F)) %>%
  'colnames<-'(c("Sample_Name", "pvalue"))
# subset beta
beta_sub <- HR_gene_list %>%
                   map_df(~(data.frame(gene = .x,
                                       CpG_probe_by_gene(., genecodefile)))) %>%
                   left_join(beta_matrix, by="CpG") %>%
                   # create columns with CpG name and corresponding values
                   pivot_longer(cols = -c("CpG", "chr_start", "chr_end", "gene", "probe_strand", "distToTSS", "CpG_context"), names_to = "Sample_Name", values_to = "beta") %>%
                  mutate(Sample_Name = if_else(startsWith(Sample_Name, "GSM"),substr(Sample_Name, 1, 10) , Sample_Name)) %>%
                   # add sample info
                   left_join(sample_sheet, by="Sample_Name")

# annotate subset beta file with CpG genome context around HR gene
gene_promo <- read.table(gene_promoter_info, header= T)

##add a column with pre-defined CpG, these probes were manually assigned from variable assessment
beta_sub <-beta_sub %>%
  mutate(PredefinedCpG = case_when(
    chr_start >= 43125041 & chr_end <= 43125714 ~ "BRCA1",
    chr_start >= 58692529 & chr_end <= 58693064 ~ "RAD51C"
  ))

##file only contains putative promoter starting, and ending location
gene_promo <- gene_promo %>%
  mutate(EPD_promoter = rownames(gene_promo),
         Start_200 = Start - 200,
         End_200 = End + 200)


##setting data.table for unequal joining
setDT(beta_sub)
setDT(gene_promo)

##unequal joining for gene annotation
beta_sub[gene_promo, on =.(chr_start > Start, chr_end < End), EPD_promo := EPD_promoter]
beta_sub[gene_promo, on =.(chr_start > Start_200, chr_end < End_200), gene_promo_200 := Gene]


#save beta file to the processed data
fwrite(as.data.table(beta_sub, keep.rownames = ""), paste0(directory, "/data/processed/methylation/",study_name, "_beta_subset_promo.csv"), row.names = F, sep = "\t")

#deduplicate samples
deduplicate_sample <- sample_sheet %>%
    mutate(Patient.ID = if_else(Project == "ICGC_OV", substr(x =  Sample.ID, start = 1, stop = 8),
                                substr(x =  Sample.ID, start = 1, stop = 15))) %>%
    dplyr::select(Sample_Name, Sample_Group, Project, Sample.ID, Patient.ID) %>%
    left_join(pvalue, by = "Sample_Name") %>%
    group_by(Patient.ID, Sample_Group, Project) %>%
    filter(pvalue == max(pvalue)) %>%
    ungroup() %>%
    dplyr::select(Sample.ID, Patient.ID, Sample_Group)

fwrite(as.data.table(deduplicate_sample, keep.rownames = ""), paste0(directory, "/data/processed/methylation/",study_name, "_deduplicated_sample.csv"), row.names = F, sep = "\t")

# beta line plot
beta_plot <- beta_plot_context_annotat(beta_sub, study_name)
beta_plot_saver(beta_plot, group = length(unique(beta_sub$Sample_Group)), output_dir = paste0(directory,"/plots/methylation/HR_gene_line"))


# beta duplicated sample plot
duplicated_sample_plot <- sample_sheet %>%
  mutate(Patient.ID = if_else(Project == "ICGC_OV",
                              substr(x =  Sample.ID, start = 1, stop = 8),
                              substr(x =  Sample.ID, start = 1, stop = 15))) %>%
  group_by(Patient.ID, Sample_Group) %>%
  filter(n() >1) %>%
  ungroup() %>%
  dplyr::select(Sample_Name, Sample.ID, Patient.ID) %>%
  left_join(beta_sub) %>%
  mutate(CpG = factor(CpG, levels = unique(beta_sub$CpG[order(beta_sub$chr_start)]))) %>%
  ggplot( aes(x = CpG, y = beta, group = Sample_Name, colour = Sample_Name)) +
              geom_point(alpha = 0.7) +
              geom_line(alpha = 0.3) +
              facet_grid(rows = vars(Patient.ID), cols = vars(gene), space = "free", scales = "free") +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5),legend.position = "none") +
              labs(title =  paste(study_name, "duplicated comparison", sep=" "))
ggsave(plot = duplicated_sample_plot,filename = paste0(directory, "/plots/methylation/HR_gene_line/", study_name,"_duplicated_sample.pdf"), height = 40, width = 30)
