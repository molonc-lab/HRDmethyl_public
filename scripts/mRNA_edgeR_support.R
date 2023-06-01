# predict whether the sample is methylated at BRCA1/RAD51C based on our criteria
BRCA1_RAD51C_meth_annotat <- function(beta_sub) {
  methylated <- beta_sub %>%
    filter(PredefinedCpG %in% c("RAD51C", "BRCA1")) %>%
    group_by(gene, Sample_Name) %>%
    summarise(probe_num = n(), brca_methylated_num = sum(beta >0.25), rad51c_methylated_num = sum(beta >0.25), fraction = brca_methylated_num/probe_num) %>%
    filter( (gene == "BRCA1" & fraction >0.6)| (gene =="RAD51C" & rad51c_methylated_num > 2)) %>%
    ungroup() %>%
    dplyr::select(gene, Sample_Name)
  beta_sub <- beta_sub %>% 
    distinct() %>%
    filter(gene %in% c("BRCA1", "RAD51C")) %>%
    mutate(BRCA1_M = ifelse(Sample_Name %in% 
                              unlist(methylated %>%filter(gene =="BRCA1") %>%dplyr::select(Sample_Name)), "Yes", "No"),
           RAD51C_M= ifelse(Sample_Name %in% 
                              unlist(methylated %>%filter(gene =="RAD51C") %>%dplyr::select(Sample_Name)), "Yes", "No"))
  return(beta_sub)
}
