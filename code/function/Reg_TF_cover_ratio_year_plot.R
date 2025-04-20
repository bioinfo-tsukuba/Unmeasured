Reg_TF_cover_ratio_year_plot <- function(calc_opt){
  
  library(tidyverse)
  library(ggrepel)
  library(RColorBrewer)
  
  tgt_region <- 500
  tgt_threshold <- 1
  
  ct_list <- read_tsv("/Users/saeko/Unmeasured/data/cell_type_list_RNA_ChIPAtlas.txt", col_names = F)
  colnames(ct_list) <- c("ct", "ct2")
  
  if(calc_opt == TRUE){
    df_binded1 <- c()
    for (tgt_row in 1:nrow(ct_list)) {
      tgt_ct <- ct_list[tgt_row, "ct"] %>% as.character() #RNA-seq
      tgt_ct2 <- ct_list[tgt_row, "ct2"] %>% as.character() #ChIP-Atlas
      
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
    write_tsv(df_binded1, "/Users/saeko/Unmeasured/data/Reg_TF_cover_ratio/Reg_TF_cover_ratio_15ct_binded.tsv")
    
  }else{
    df_binded1 <- read_tsv("/Users/saeko/Unmeasured/data/Reg_TF_cover_ratio/Reg_TF_cover_ratio_15ct_binded.tsv")
    df_binded2 <- df_binded1 %>% mutate(label = ifelse(year == 2011 | year == 2020, Cell_type, ""))
    color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG"))
    
    p1 <- df_binded2 %>% ggplot(aes(x = year, y = RR_cover_ratio, color = Cell_type, label = label)) +
      geom_point() +
      geom_line() +
      geom_label_repel(size = 10) +
      scale_color_manual(values = color_list) +
      xlab("Year")+
      ylab("Reg-TF cover ratio")+
      theme(plot.title = element_text(size = unit(15, "pt"), face="bold"), 
            #legend.position = c(0.12,0.73),
            legend.position = "none",
            legend.title = element_text(size = unit(10, "pt")),
            legend.text = element_text(size = unit(10, "pt")),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.ticks = element_line(colour = "black", linewidth = unit(0.5, "pt")),
            axis.line = element_line(linewidth = unit(0.5, "pt")),
            axis.title = element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
            axis.text = element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
            axis.text.x =element_text(size = unit(30, "pt"), colour = "black", face = "bold", angle = 45, hjust = 1),
            # axis.line = element_line(colour="black"),
            aspect.ratio = 1
      ) 
      
    
    df_binded3 <- df_binded1 %>% mutate(label = ifelse(year == 2011, Cell_type, ""))
    p2 <- df_binded3 %>% 
      ggplot(aes(x = log(uniq_TF_num, base = 10), y = RR_cover_ratio, color = Cell_type, label = label)) +
      geom_point() +
      geom_line()+
      geom_label_repel(size = 10) +
      scale_color_manual(values = color_list) +
      xlab("Number of unique TFs (log10)")+
      ylab("Reg-TF cover ratio")+
      theme(plot.title = element_text(size = unit(15, "pt"), face="bold"), 
            legend.position = c(0.82,0.45),
            #legend.position = "none",
            legend.title = element_blank(),
            legend.text = element_text(size = unit(20, "pt")),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.ticks = element_line(colour = "black", linewidth = unit(0.5, "pt")),
            axis.line = element_line(linewidth = unit(0.5, "pt")),
            axis.title = element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
            axis.text = element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
            axis.text.x =element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
            # axis.line = element_line(colour="black"),
            aspect.ratio = 1
      ) 
  }
  
  return(list(p1, p2))
}