library(methylKit)
library(tidyverse)
library(ggplot2)
library(dplyr)

mapping <- c(
  "SRR8633864"="MKN74_STOMACH",
  "SRR8633865"="MKN7_STOMACH",
  "SRR8633883"="SNU668_STOMACH",
  "SRR8633887"="SNU601_STOMACH",
  "SRR8633889"="SNU520_STOMACH",
  "SRR8633890"="SNU5_STOMACH",
  "SRR8633902"="MKN1_STOMACH",
  "SRR8633908"="MKN45_STOMACH",
  "SRR8633910"="IM95_STOMACH",
  "SRR8633941"="HGC27_STOMACH",
  "SRR8633967"="SNU719_STOMACH",
  "SRR8633973"="SNU620_STOMACH",
  "SRR8634018"="RERFGC1B_STOMACH",
  "SRR8634035"="KE39_STOMACH",
  "SRR8634072"="GSS_STOMACH",
  "SRR8634073"="GSU_STOMACH",
  "SRR8634111"="NUGC4_STOMACH",
  "SRR8633575"="SNU16_STOMACH",
  "SRR8633579"="NUGC3_STOMACH",
  "SRR8633591"="TGBC11TKB_STOMACH",
  "SRR8633615"="FU97_STOMACH",
  "SRR8633634"="2313287_STOMACH",
  "SRR8633673"="SH10TC_STOMACH",
  "SRR8633707"="KATOIII_STOMACH",
  "SRR8633738"="GCIY_STOMACH",
  "SRR8633801"="SNU1_STOMACH",
  "SRR8633802"="SNU216_STOMACH",
  "SRR8633850"="NUGC2_STOMACH",
  "SRR8633212"="HUG1N_STOMACH",
  "SRR8633294"="AGS_STOMACH",
  "SRR8633310"="LMSU_STOMACH",
  "SRR8633365"="ECC10_STOMACH",
  "SRR8633366"="ECC12_STOMACH",
  "SRR8633424"="NCIN87_STOMACH",
  "SRR8633455"="HS746T_STOMACH",
  "SRR8633471"="NCCSTCK140_STOMACH",
  "SRR8633526"="OCUM1_STOMACH"
)


file_list = paste0("/working/lab_olgak/brettL/work/stad_depmap/rerun/data/intermediate_data/bis_aligned_sorted/", names(mapping),"_bismark_bt2.sorted.bam")

treat_vec <- rep(1, length(file_list))

list <- file_list[1:37] %>%
  map(~ processBismarkAln(.x , assembly = "hg37", sample.id= gsub("_bismark_bt2.sorted.bam", "",basename(.x)),mincov = 5, minqual = 20))

save.image("/working/lab_olgak/brettL/work/stad_depmap/rerun/analysis/rdata_files/mincov5minqual20STADDEPMAP.RData")
