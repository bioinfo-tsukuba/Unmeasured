TFs_per_fam <- function(){
  
  library(tidyverse)
  # 20231004 version
  df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/SRX_date_20231004.tsv")
  df <- df %>% arrange(year, month, day, time) %>% mutate(order = 1:nrow(df)) %>% 
    filter(Cell_type_class != "Unclassified" & Cell_type_class != "No description") %>%
    filter(Antigen != "GFP" & Antigen != "Epitope tags")
  
  df_fam <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS_paper_public/data/Fig2/obj/df_fam.rds")
  tgt_IDs <- df$SRX %>% unique() %>% as.character()
  tgt_df_fam <- df_fam %>% filter(ID %in% tgt_IDs) %>% select(Antigen, Family) %>% distinct()
  
  df3 <- df %>% left_join(tgt_df_fam, by = "Antigen") %>% 
    filter(is.na(Family) == FALSE & Family != "No_annotation")
  
  p <- df3 %>% group_by(Family) %>% summarise(uniq_TF = n_distinct(Antigen)) %>%
    ggplot(aes(x = reorder(Family, -uniq_TF), y = uniq_TF)) +
    geom_bar(stat = "identity") +
    xlab("TF Family")+
    ylab("Total number of unique TFs")+
    #scale_color_manual(values = color_list) +
    theme(plot.title = element_text(size = unit(15, "pt"), face="bold"), 
          legend.position = "none",
          legend.title = element_text(size = unit(15, "pt")),
          legend.text = element_text(size = unit(15, "pt")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks = element_line(colour = "black", linewidth = unit(0.5, "pt")),
          axis.line = element_line(linewidth = unit(0.5, "pt")),
          axis.title = element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
          axis.text = element_text(size = unit(30, "pt"), colour = "black", face = "bold"),
          axis.text.x =element_text(size = unit(10, "pt"), colour = "black", face = "bold", angle = 45, hjust = 1),
          # axis.line = element_line(colour="black"),
          aspect.ratio = 1
    )
  return(p)
}