# predict whether the sample is methylated at BRCA1/RAD51C based on our criteria
BRCA1_RAD51C_meth_annotat <- function(beta_sub) {
  methylated <- beta_sub %>%
    filter(PredefinedCpG %in% c("RAD51C", "BRCA1")) %>%
    group_by(gene, Sample_Name) %>%
    summarise(probe_num = round(n()* 0.6, 0), methylated_num = sum(round(beta,2) >= 0.25)) %>%
    filter(probe_num <= methylated_num) %>%
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
