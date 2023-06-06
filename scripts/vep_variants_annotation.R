args = commandArgs(trailingOnly=TRUE)

working_dir = args[1]

# TCGA in a single file MC3 currated
tcga_mc3_location = paste0(working_dir, "/data/raw/mc3.v0.2.8.PUBLIC.maf")
# ICGC in separate MAF files we use somatic high confidence only, here a folder contain all the MAF files will be suffice
icgc_mc3_location = paste0(working_dir, "/data/raw/20220505_ICGC_OV_meth/MAFs/")


out_dir = paste0(working_dir, "/data/processed/VCF")

library("data.table")
library("tidyverse")

tcga_mc3 <- as.data.frame(fread(tcga_mc3_location))
icgc_mc3 <- Sys.glob(paste0(icgc_mc3_location, "/*shc.maf.gz")) %>%
  map_dfr(~fread(cmd = paste0("gunzip -cq ", .x), skip = "Hugo_Symbol") %>%
            filter(Hugo_Symbol != "unknown"))

hr_gene = c("TP53", "BRCA1", "BRCA2", "CDK12", "RAD51C", "PALB2", "BRIP1", "RAD51D")
gene_list = c("TP53", "PIK3CA", "IDH1", "ARID1A", "PTEN", "NUMA1", "QKI", "BCL11B", "BAP1", "CIC", "ATRX")

cat("##fileformat=VCFv4.1\n##contig=<ID=1,length=249250621,assembly=b37>\n##reference=file:///path/to/human_g1k_v37.fasta\n",
    append=FALSE,
    file= paste0(out_dir, "/variants_hg19.vcf"))

HR_tp53_vcf <- rbind(tcga_mc3 %>%
                       filter(Hugo_Symbol %in% unique(c(hr_gene, gene_list))) %>%
                       mutate(`#CHROM` = str_remove(Chromosome,"chr"),
                              POS = Start_Position,
                              ID = ".",
                              REF = Reference_Allele,
                              ALT = if_else(Tumor_Seq_Allele1 == Reference_Allele, Tumor_Seq_Allele2, Tumor_Seq_Allele1 ),
                              QUAL = ".",
                              FILTER = "PASS",
                              INFO = ".") %>%
                       select(`#CHROM`,POS,ID,REF,ALT,QUAL,FILTER,INFO),
                     icgc_mc3%>%
                       filter(Hugo_Symbol %in% unique(c("TP53", "BRCA1", "BRCA2", "CDK12", "RAD51C", "PALB2", "BRIP1", "RAD51D", gene_list)))%>%
                       mutate(`#CHROM` = Chromosome,
                              POS = Start_Position,
                              ID = ".",
                              REF = Reference_Allele,
                              ALT = if_else(Tumor_Seq_Allele1 == Reference_Allele, Tumor_Seq_Allele2, Tumor_Seq_Allele1 ),
                              QUAL = ".",
                              FILTER = "PASS",
                              INFO = ".") %>%
                       select(`#CHROM`,POS,ID,REF,ALT,QUAL,FILTER,INFO)) %>%
  arrange(`#CHROM`, POS) %>%
  distinct()

fwrite(HR_tp53_vcf,
       append=TRUE,
       sep="\t",
       col.names = T,
       file=paste0(out_dir, "/variants_hg19.vcf"))
