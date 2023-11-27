UniqueTF_RRratio_scatter <- function(tgt_region){
  
  library(tidyverse)
  library(RColorBrewer)
  library(patchwork)
  
  ct_list <- read_tsv("/Users/saeko/Unmeasured/data/cell_type_list_RNA_ChIPAtlas.txt", col_names = F)
  colnames(ct_list) <- c("ct", "ct2")
  
  patch <- c()
  for (tgt_threshold in 1:10) {
    
    print(tgt_threshold)
    df_binded1 <- c()
    for (tgt_row in 1:nrow(ct_list)) {
      tgt_ct <- ct_list[tgt_row, "ct"] %>% as.character() #RNA-seq
      tgt_ct2 <- ct_list[tgt_row, "ct2"] %>% as.character() #ChIP-Atlas
      tgt_region <- 500
      
      print(tgt_ct)
      print(tgt_ct2)
      
      RR_cover_ratio_tib <- suppressMessages(read_tsv(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/RR_ratio_table/RR_cover_ratio_", tgt_ct2, "_",tgt_region, ".tsv")))
      uniq_TF_df <- suppressMessages(read_tsv(paste0("/Users/saeko/Unmeasured/data/uniqTF_table/", tgt_ct2, "_", tgt_region, "bp.tsv")))
      tgt_RR_cover_ratio_tib <- RR_cover_ratio_tib %>% filter(threshold == tgt_threshold)
      df_join1 <- tgt_RR_cover_ratio_tib %>% left_join(uniq_TF_df, by = "year") %>%
        mutate(Cell_type = tgt_ct2, tgt_region = tgt_region)
      
      if(tgt_row == 1){
        df_binded1 <- df_join1
      }else{
        df_binded1 <- rbind(df_binded1, df_join1)
      }
    } #for, tgt_row
    
    color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(5,"BrBG"))
    p1 <- df_binded1 %>% ggplot(aes(x = uniq_TF_num, y = RR_cover_ratio, color = Cell_type)) +
      geom_point() +
      geom_line()+
      scale_color_manual(values = color_list) +
      ggtitle(paste0("Threshold = ", tgt_threshold)) +
      xlab("Number of unique TFs")+
      ylab("Ratio")+
      theme(plot.title = element_text(face="bold",hjust = 0.5), 
            legend.position = "none",
            #legend.position = c(0.9, 0.4),
            legend.background = element_rect(color = "black", fill = alpha("white", alpha = 1)),
            panel.grid.major = element_line(colour="gray"),
            panel.grid.minor = element_line(colour="gray", size = 1),
            panel.background = element_blank(), 
            axis.line = element_line(colour="black"),
            axis.text=element_text(size=3,face="bold", color = "black"),
            axis.text.x =element_text(size=3,face="bold", color = "black"),
            axis.text.y =element_text(size=3,face="bold", color = "black"),
            axis.title=element_text(size=3,face="bold", color = "black"),
            aspect.ratio = 0.5
      )
    
    p2 <- df_binded1 %>% ggplot(aes(x = uniq_TF_num, y = RR_cover_ratio, color = Cell_type)) +
      geom_point() +
      geom_line()+
      scale_color_manual(values = color_list) +
      ggtitle(paste0("Threshold = ", tgt_threshold)) +
      xlab("Number of unique TFs")+
      ylab("Ratio")+
      theme(plot.title = element_text(face="bold",hjust = 0.5), 
            #legend.position = "none",
            #legend.position = c(0.9, 0.4),
            legend.background = element_rect(color = "black", fill = alpha("white", alpha = 1)),
            panel.grid.major = element_line(colour="gray"),
            panel.grid.minor = element_line(colour="gray", size = 1),
            panel.background = element_blank(), 
            axis.line = element_line(colour="black"),
            axis.text=element_text(size=15,face="bold", color = "black"),
            axis.text.x =element_text(size=15,face="bold", color = "black"),
            axis.text.y =element_text(size=15,face="bold", color = "black"),
            axis.title=element_text(size=15,face="bold", color = "black"),
            aspect.ratio = 0.5
      )
    ggsave(paste0("/Users/saeko/Unmeasured/plot/UniqueTF_RRratio_scatter/tmp/UniqueTF_RRratio_scatter_", tgt_ct2, "_threshold_", tgt_threshold, ".pdf"), p2, width = 9, height = 7)
    
    if(tgt_threshold == 1){
      patch <- p1
    }else{
      patch <- patch + p1
    }
  }
  return(list(patch))
}