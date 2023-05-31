idat_name_generator <- function(filename, type){
  filname <- gsub("_Grn.idat|_Red.idat", "", filename)
  name_list = unlist(str_split(filname, "_"))
  if( length(name_list) >2) {
    id = paste(name_list[1:length(name_list) - 1],collapse = '_')
    position = name_list[length(name_list)]
  } else{
    id = name_list[1]
    position = name_list[2]
  }
  if (toupper(type) == "ID") {
    return(id)
  } else if (toupper(type) == "POSITION") {
    return(position)
  }
}

#read geo matrix to get sample name mapping
load_sample_geo <- function(filepath, keywords){
  #read sample types
  line <- fread(filepath, header = F, sep = "\t", skip = keywords, nrows = 1 )
  Sample_Name <- fread(filepath, header = F, sep = "\t", skip = "ID_REF", nrows = 1 )
  
  #fix formatting and typo
  line <- gsub("tumour","tumor", line, fixed = T)
  line <- gsub("primarytumour","primarytumor", line, fixed = T)
  line <- tolower(line)
  # check for the line read to be containing sample type
  if(TRUE %in% str_detect(line, "tumor|normal|primary|n$|t$|N$|T$")){
    sample_type <- as.data.frame(cbind(unlist(Sample_Name[,-1]),unlist(line[2:length(line)])))
    colnames(sample_type)<-c("Sample_Name","Sample_Group")
  }
  
  sample_type$Sample_Group <- gsub("umor", "umour", sample_type$Sample_Group)
  #summarise sample type
  print(sample_type%>%
          dplyr::group_by(Sample_Group) %>%
          dplyr::summarise(n = n())
  )
  return(sample_type)
}

# generate pd file for champ to run IDAT
# input sample file: e.g brca_sample.tsv, this file could be obtained from TCGA when getting manifest file
#
tcga_pd_generator <- function(sample_file) {
  # load sample file
  study_name = paste(unlist(strsplit(sample_file, split = "_"))[2:3], collapse ="_")
  sample_pd <- fread(sample_file, header = T, sep = '\t')
  out_dir = dirname(sample_file)
  
  sample_pd <- sample_pd %>% 
    mutate(Sample_Name = unlist(lapply(`File Name`, idat_name_generator, "id")),
           Sample_Plate = NA, 
           Sample_Group = `Sample Type`, 
           Pool_ID = NA, 
           Project = study_name, 
           Sample_Well = NA, 
           Sentrix_ID = unlist(lapply(`File Name`, idat_name_generator,"id")), Sentrix_Position = unlist(lapply(`File Name`, idat_name_generator,"position"))) %>%
    dplyr::select(Sample_Name,Sample_Plate, Sample_Group, Pool_ID, Project, Sample_Well, Sentrix_ID, Sentrix_Position, Sample_Group, Sample.ID = `Sample ID`) %>%
    distinct()
  
  write.csv(sample_pd, paste0(out_dir, "/IDAT/", study_name, "_pd.csv"), row.names = FALSE)
}

geo_pd_generator <- function(sample_matrix, keywords) {
  study_name = paste(unlist(strsplit(sample_matrix, split = "_"))[2:3], collapse ="_")
  out_dir = dirname(sample_matrix)
  
  sample_pd <- Sys.glob(paste0(dirname(sample_matrix), "/IDAT/*.idat"))%>%
    as.data.frame() %>%
    mutate(File = basename(.)) %>%
    separate(File, into = c("V1", "V2", "V3"), extra = "drop") %>%
    dplyr::select("V1", "V2", "V3")
  
  sample_type = load_sample_geo(sample_matrix,keywords)
  
  sample_pd <- left_join(sample_pd, sample_type, by = c("V1" = "Sample_Name")) %>%
    mutate(Sample_Name = V1, Sample_Plate = NA, Pool_ID = NA, Project = study_name, Sample_Well = NA, Sentrix_ID = V2, Sentrix_Position= V3) %>%
    dplyr::select(Sample_Name,Sample_Plate, Sample_Group, Pool_ID, Project, Sample_Well, Sentrix_ID, Sentrix_Position) %>%
    distinct()
  return(sample_pd)
}

