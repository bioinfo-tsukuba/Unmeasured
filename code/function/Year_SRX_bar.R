Year_SRX_bar <- function(){
  
  library(tidyverse)
  library(ggrepel)
  library(RColorBrewer)
  
  # 20231004 version
  df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/SRX_date_20231004.tsv")
  df <- df %>% arrange(year, month, day, time) %>% mutate(order = 1:nrow(df)) %>% 
    filter(Cell_type_class != "Unclassified" & Cell_type_class != "No description") %>%
    filter(Antigen != "GFP" & Antigen != "Epitope tags")
  
  color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG"))
  p1 <- df %>% 
    #group_by(year) %>% 
    #summarise(SRX_num =n()) %>%
    ggplot(aes(x = factor(year), fill = Cell_type_class))+
    geom_bar() +
    scale_fill_manual(values = color_list) +
    xlab("Year")+
    ylab("Total number of SRX (TF ChIP-seq)")+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          #legend.position = "none",
          panel.grid.major = element_line(colour="gray"),
          panel.grid.minor = element_line(colour="gray", size = 1),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=15,face="bold", color = "black"),
          axis.text.x =element_text(size=15,face="bold", color = "black", angle = 45, hjust = 1),
          axis.text.y =element_text(size=15,face="bold", color = "black"),
          axis.title=element_text(size=15,face="bold", color = "black"),
          aspect.ratio = 0.7
    )
  
}