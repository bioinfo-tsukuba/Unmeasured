function(tgt_region){
  
  Low_map_genes_df <- read_tsv("/Users/saeko/Unmeasured/data/ENCODE_blacklist/Low_mappability_genes_geneID.txt", col_names = F)
  Low_map_genes <- Low_map_genes_df$X2 %>% unique()
  
  ct_list <- read_tsv("/Users/saeko/Unmeasured/data/cell_type_list_RNA_ChIPAtlas.txt", col_names = F)
  colnames(ct_list) <- c("ct", "ct2")
  patch1 <- c()
  Reg_TF_ratio_ceiling_list <- list()
  for (tgt_row in 1:nrow(ct_list)) {
    tgt_ct <- ct_list[tgt_row, "ct"] %>% as.character() #RNA-seq
    tgt_ct2 <- ct_list[tgt_row, "ct2"] %>% as.character() #ChIP-Atlas
    df_wider_selected_annotated <- suppressMessages(read_tsv(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/RR_ratio_table/df_wider_selected_annotated_", tgt_ct2, "_",tgt_region, "bp.tsv")))
    RR_cover_ratio_tib <- suppressMessages(read_tsv(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/RR_ratio_table/RR_cover_ratio_", tgt_ct2, "_",tgt_region, ".tsv")))
    
    all_expressed_gene <- colnames(df_wider_selected_annotated[,12:ncol(df_wider_selected_annotated)])
    Low_map_expressed_genes <- intersect(all_expressed_gene, Low_map_genes)
    Reg_TF_ratio_ceiling <- 1 - length(Low_map_expressed_genes) / length(all_expressed_gene)
    print(Low_map_expressed_genes)
    print(paste0(tgt_ct2, " : ", Reg_TF_ratio_ceiling))
    Reg_TF_ratio_ceiling_list[[tgt_ct2]] <- Reg_TF_ratio_ceiling
    
    color_list <- c(brewer.pal(10,"Spectral"))
    p1 <- RR_cover_ratio_tib %>% 
      filter(threshold %in% 1:10) %>%
      ggplot(aes(x = year, y = RR_cover_ratio, color = factor(threshold))) +
      geom_point() +
      geom_line() +
      ylim(c(0,1))+
      geom_hline(yintercept=Reg_TF_ratio_ceiling, linetype=2,alpha=0.7,size=1,color='red4') +
      xlab("Year")+
      ylab("Ratio")+
      labs(color = "Threshold")+
      scale_color_manual(values = color_list) +
      #ggtitle(paste0("Regulatory-TF cover ratio in expressed genes, ", tgt_ct2))+
      ggtitle(tgt_ct2) +
      theme(plot.title = element_text(face="bold",hjust = 0.5), 
            legend.position = "none",
            #legend.position = c(0.9, 0.4),
            legend.background = element_rect(color = "black", fill = alpha("white", alpha = 1)),
            panel.grid.major = element_line(colour="gray"),
            panel.grid.minor = element_line(colour="gray", size = 1),
            panel.background = element_blank(), 
            axis.line = element_line(colour="black"),
            axis.text=element_text(size=3,face="bold", color = "black"),
            axis.text.x =element_text(size=3,face="bold", color = "black", angle = 45, hjust = 1),
            axis.text.y =element_text(size=3,face="bold", color = "black"),
            axis.title=element_text(size=3,face="bold", color = "black"),
            aspect.ratio = 0.5
      )
    ggsave(paste0("/Users/saeko/Unmeasured/plot/RR_ratio_year_ceiling/", tgt_ct2, "_", tgt_region, "bp_ceiling.pdf" ), p1, width = 9, height = 7)
    
    if(tgt_row == 1){
      patch1 <- p1
    }else{
      patch1 <- patch1 + p1
    }
  }
  
  ggsave(paste0("/Users/saeko/Unmeasured/plot/RR_ratio_year_ceiling/allTF_", tgt_region, "bp.pdf" ), patch1)
  
  Reg_TF_ratio_ceiling_df <- tibble(Cell_type = names(Reg_TF_ratio_ceiling_list), ceiling = as.numeric(Reg_TF_ratio_ceiling_list))
  Reg_TF_ratio_ceiling_df %>% ggplot(aes(x = reorder(Cell_type, -ceiling), y = ceiling, label = ceiling)) +
    geom_point()+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = "none",
          #legend.position = c(0.9, 0.4),
          legend.background = element_rect(color = "black", fill = alpha("white", alpha = 1)),
          panel.grid.major = element_line(colour="gray"),
          panel.grid.minor = element_line(colour="gray", size = 1),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=10,face="bold", color = "black"),
          axis.text.x =element_text(size=10,face="bold", color = "black", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold", color = "black"),
          axis.title=element_text(size=10,face="bold", color = "black"),
          aspect.ratio = 0.5
    )
  
  
  
}