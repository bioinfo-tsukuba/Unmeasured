Uniq_TF_year <- function(tgt_ct, tgt_ct2, tgt_region){
  
  library(tidyverse)
  library(RColorBrewer)
  print(tgt_ct)
  print(tgt_ct2)
  
  df_date <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/SRX_date_20231004.tsv")
  df_date <- df_date %>% arrange(year, month, day, time) %>% mutate(order = 1:nrow(df_date)) %>% 
    filter(Cell_type_class != "Unclassified" & Cell_type_class != "No description") %>%
    filter(Antigen != "GFP" & Antigen != "Epitope tags")
  df_date_selected <- df_date %>% filter(Cell_type == tgt_ct2)
  
  
  # Uniq TF - year -----
  year_list_input <- unique(df_date$year) %>% sort()
  uniq_TF_list <- c()
  year_list <- c()
  for(tgt_year in year_list_input) {
    print(tgt_year)
    df_date2 <- df_date_selected %>% filter(year %in% year_list_input[1]:tgt_year)
    uniq_TF_list <- c(uniq_TF_list, length(unique(df_date2$Antigen)))
    year_list <- c(year_list, tgt_year)
  }
  uniq_TF_df <- tibble(year = year_list, uniq_TF_num = uniq_TF_list)
  
  p1 <- uniq_TF_df %>% ggplot(aes(x = year, y = uniq_TF_num)) +
    geom_point()+
    geom_line()+
    xlab("Year")+
    ylab("Number of unique TFs")+
    ggtitle(paste0("Number of unique TFs in ", tgt_ct2))+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          #legend.position = "none",
          #legend.position = c(0.9, 0.4),
          legend.background = element_rect(color = "black", fill = alpha("white", alpha = 1)),
          panel.grid.major = element_line(colour="gray"),
          panel.grid.minor = element_line(colour="gray", size = 1),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=15,face="bold", color = "black"),
          axis.text.x =element_text(size=15,face="bold", color = "black", angle = 45, hjust = 1),
          axis.text.y =element_text(size=15,face="bold", color = "black"),
          axis.title=element_text(size=15,face="bold", color = "black"),
          aspect.ratio = 0.5
    )
  write_tsv(uniq_TF_df, paste0("/Users/saeko/Unmeasured/data/uniqTF_table/", tgt_ct2, "_", tgt_region, "bp.tsv"))
  
  # SRX number - year ------
  year_list_input <- unique(df_date$year) %>% sort()
  SRX_num_list <- c()
  year_list <- c()
  for(tgt_year in year_list_input) {
    print(tgt_year)
    df_date2 <- df_date_selected %>% filter(year %in% year_list_input[1]:tgt_year)
    SRX_num_list <- c(SRX_num_list, length(unique(df_date2$SRX)))
    year_list <- c(year_list, tgt_year)
  }
  SRX_df <- tibble(year = year_list, SRX_num = SRX_num_list)
  
  p2 <- SRX_df %>% ggplot(aes(x = year, y = SRX_num)) +
    geom_point()+
    geom_line()+
    xlab("Year")+
    ylab("Number of ChIP-seq")+
    ggtitle(paste0("Number of ChIP-seq in ", tgt_ct2))+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          #legend.position = "none",
          #legend.position = c(0.9, 0.4),
          legend.background = element_rect(color = "black", fill = alpha("white", alpha = 1)),
          panel.grid.major = element_line(colour="gray"),
          panel.grid.minor = element_line(colour="gray", size = 1),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=15,face="bold", color = "black"),
          axis.text.x =element_text(size=15,face="bold", color = "black", angle = 45, hjust = 1),
          axis.text.y =element_text(size=15,face="bold", color = "black"),
          axis.title=element_text(size=15,face="bold", color = "black"),
          aspect.ratio = 0.5
    )
  write_tsv(SRX_df, paste0("/Users/saeko/Unmeasured/data/SRX_table/", tgt_ct2, "_", tgt_region, "bp.tsv"))
  
  return(list(p1, p2))
  
  
}