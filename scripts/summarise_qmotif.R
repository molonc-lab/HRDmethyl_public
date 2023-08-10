library(tidyverse)
library(xml2)
library(readxl)
metasheet="../data/processed/AOCS_-_ICGC_ovarian_cancer_project_qMotif_Report.xlsx"


# List all qmotif output files
filenames_qmotif <- data.frame(read_xlsx(metasheet)) %>%
  dplyr::select(DonorID = `Donor.Label`, SampleID = `CollectedSample.Label`, file_path = `qMotif.Analysis.FilePath`)
filenames_qmotif <- list.files(path=paste0(path,"output"),
                        full.names = T,
                        pattern = ".xml$",
                        recursive = F)
  
all_counts <- filenames_qmotif %>%
  pull(file_path)%>%
  map_dfr(~ read_xml(Sys.glob(paste0(.x, "/*.xml"))) %>% 
            xml_find_all("//summary/counts") %>% 
            xml_children() %>%
            map_df(~list(count=xml_attr(.,"count"),group=xml_name(.))) %>%
            spread(group,count) %>%
            mutate(DonorID = filenames_qmotif[which(filenames_qmotif$file_path == .x),"DonorID"],
                   SampleID = filenames_qmotif[which(filenames_qmotif$file_path == .x),"SampleID"]))


write.table(all_counts, paste0("../data/processed/","qmotif_all_counts.tsv"), sep = "\t", col.names = T, quote = F, row.names = F)





