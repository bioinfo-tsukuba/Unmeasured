Ct_per_Ctc <- function(){
  
  library(tidyverse)
  # 20231004 version
  df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/SRX_date_20231004.tsv")
  df <- df %>% arrange(year, month, day, time) %>% mutate(order = 1:nrow(df)) %>% 
    filter(Cell_type_class != "Unclassified" & Cell_type_class != "No description") %>%
    filter(Antigen != "GFP" & Antigen != "Epitope tags")
  
  p <- df %>% group_by(Cell_type_class) %>% summarise(uniq_ct = n_distinct(Cell_type)) %>%
    ggplot(aes(x = reorder(Cell_type_class, -uniq_ct), y = uniq_ct)) +
    geom_bar(stat = "identity") +
    xlab("Cell type class")+
    ylab("Total number of unique Cell type")+
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