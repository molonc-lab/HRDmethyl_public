CpG_probe_by_gene<-function(gene_name, genecodefile){
  #default download gz file
  if (grepl(".gz", genecodefile)){
    gene_table <- fread(cmd = paste("gunzip -cq ", genecodefile, " |", " grep ", toupper(gene_name), sep = ""), header = F, sep = "\t")
  } else {
    # allow to feed in txt, csv, tsv file if needed
    gene_table <- fread(cmd = paste("grep",toupper(gene_name), genecodefile, sep = " "),sep = "\t")
  }
  gene_table <- gene_table %>%
    dplyr::select("V2", "V3", "V4", "V5", "V11", "V12") %>%
    dplyr::rename(chr_start = V2, chr_end = V3, probe_strand = V4, CpG = V5, distToTSS = V11, CpG_context = V12)
  return(gene_table)
}

#plot beta line plot
beta_plot_context_annotat <- function(beta_sub, study_name){
  beta_plot <- beta_sub %>%
    dplyr::group_by(gene) %>%
    mutate(CpG = factor(CpG, levels = unique(beta_sub$CpG[order(beta_sub$chr_start)])), gene2=gene) %>%
    group_map(.f = ~ ggplot(.x, aes(x = CpG, y = beta, group = Sample_Name, colour = Sample_Name)) +
                geom_point(alpha = 0.7) +
                geom_line(alpha = 0.3) +
                facet_wrap(~Sample_Group, ncol = 1) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
                guides(colour = "none") +
                scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0), expand=c(0,0)) +
                
                geom_tile(aes(x= CpG, y = -0.05, fill = gene_promo_200, height = 0.05), color = NA) +
                scale_fill_brewer("Promoter +- 200bps", palette = "PRGn",  guide = guide_legend(order = 1))+
                new_scale_fill() +
                
                geom_tile(aes(x= CpG, y = -0.12, fill = CpG_context, height = 0.05), color = NA) +
                scale_fill_brewer("CpG context", palette = "Greens",  guide = guide_legend(order = 2),  direction=-1)+
                new_scale_fill() +
                geom_tile(aes(x= CpG, y = -0.19, fill = EPD_promo, height = 0.05), color = NA) +
                scale_fill_brewer("Predicted promoter", palette = "Set3",  guide = guide_legend(order = 3))+
                new_scale_fill() +
                geom_tile(aes(x= CpG, y = -0.26, fill = PredefinedCpG, height = 0.05), color = NA) +
                scale_fill_brewer("CpG selected", palette = "PiYG",  guide = guide_legend(order = 4))+
                
                
                labs(title =  paste(study_name, unique(.x$gene2), "Methylation", sep="")), .keep = T)
  names(beta_plot) <- sapply(unique(sort(beta_sub$gene)), function(x) paste(study_name, x, sep="_"))
  
  names(beta_plot) <- sapply(unique(sort(beta_sub$gene)), function(x) paste(study_name, x, sep="_"))
  return(beta_plot)
}


beta_plot_saver <- function(plots, group = 2, output_dir, plot_name = "_methylation_overview.pdf") {
  lapply(names(plots), 
         function(x) ggsave(filename = paste0(output_dir,"/", x, plot_name), plot=plots[[x]], height = 12 * group , width = 30, units = "cm"))
}
