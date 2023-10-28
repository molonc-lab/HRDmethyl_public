args = commandArgs(trailingOnly=TRUE)


library("data.table")
library("tidyverse")

working_dir = args[1]
out_dir = paste0(working_dir, "/data/processed/VCF")
maf_file = paste0(working_dir, "/data/processed/depmap_brett/OmicsSomaticMutations_23Q2_0523.csv")

gene_list <- c("BRCA1", "BRCA2", "RAD51C", "PALB2", "BRIP1", "RAD51D", "TP53", "PIK3CA", "ARID1A")

mutation_depmap <- fread(cmd = paste0("grep -E 'TP53|PIK3CA|ARID1A|BRCA|RAD51|PALB2|BRIP1' ",maf_file)) %>%
  `colnames<-`( colnames(fread(cmd = paste0("head ",maf_file)))) %>%
  filter(HugoSymbol %in% gene_list) %>%
  right_join(depmap_model %>%
               dplyr::select(ModelID, Sample.ID = StrippedCellLineName)) %>%
  filter(Sample.ID %in% depmap_mRNA$Cell_name, !VariantInfo %in% c("INTRON", "FIVE_PRIME_FLANK", "THREE_PRIME_UTR" ,"SILENT"))


#clinvar annotation
cat("##fileformat=VCFv4.1\n##contig=<ID=1,length=249250621,assembly=b38>\n##reference=file:///path/to/human_g1k_v38.fasta\n",
    append=FALSE,
        file= "data/processed/VCF/depmap_variants_hg38.vcf")

HR_tp53_vcf <- mutation_depmap %>%
                       mutate(`#CHROM` = as.numeric(gsub("chr", "", Chrom)),
                              POS = Pos,
                              ID = ".",
                              REF = Ref,
                              ALT =Alt,
                              QUAL = ".",
                              FILTER = "PASS",
                              INFO = ".") %>%
                       select(`#CHROM`,POS,ID,REF,ALT,QUAL,FILTER,INFO) %>%
  arrange(`#CHROM`, POS) %>%
  distinct()

fwrite(HR_tp53_vcf,
       append=TRUE,
       sep="\t",
       col.names = T,
       file=paste0(out_dir, "/depmap_variants_hg19.vcf"))

