comp_year <- function(){
  
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
  
  # 20231004 version
  df <- suppressMessages(read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/SRX_date_20231004.tsv"))
  df <- df %>% arrange(year, month, day, time) %>% mutate(order = 1:nrow(df)) %>% 
    filter(Cell_type_class != "Unclassified" & Cell_type_class != "No description") %>%
    filter(Antigen != "GFP" & Antigen != "Epitope tags")
  
  patch <- c()
  for (tgt_year in (min(df$year)+1):max(df$year)) {
    
    print(tgt_year)
    df1 <- df %>% filter(year %in% min(unique(df$year)):tgt_year) %>%  
      group_by(Antigen) %>% summarise(SRX_num_x = n()) %>% mutate(axis_x = "x")
    df2 <- df %>% filter(year %in% (tgt_year+1):max(unique(df$year)))%>%  
      group_by(Antigen) %>% summarise(SRX_num_y = n()) %>% mutate(axis_y = "y")
    df3 <- df1 %>% full_join(df2, by = "Antigen") 
    df3$SRX_num_x[is.na(df3$SRX_num_x)] <- 0
    df3$SRX_num_y[is.na(df3$SRX_num_y)] <- 0
    df3$axis_x[is.na(df3$axis_x)] <- "x"
    df3$axis_y[is.na(df3$axis_y)] <- "y"
    
    cor <- cor.test(log(df3$SRX_num_x, base = 10), log(df3$SRX_num_y, base = 10), method = "spearman")
    print(cor)
    
    p1 <- df3 %>% ggplot(aes(x = log(SRX_num_x, base = 10), y = log(SRX_num_y, base = 10), label = Antigen)) +
      geom_point()+
      xlab(paste0("No. of log10(ChIP-seq) until ", tgt_year)) +
      ylab(paste0("No. of log10(ChIP-seq) ", tgt_year+1, " - ", max(df$year)))+
      geom_text_repel(size = 5, force = 10) +
      annotate("text", x=log(max(df3$SRX_num_x), base = 10),   y= 0.5, label = paste0("spearman : ", cor$estimate)) +
      theme(plot.title = element_text(face="bold",hjust = 0.5), 
            #legend.position = "none",
            panel.grid.major = element_line(colour="gray"),
            panel.grid.minor = element_line(colour="gray", size = 1),
            panel.background = element_blank(), 
            axis.line = element_line(colour="black"),
            axis.text=element_text(size=15,face="bold", color = "black"),
            axis.text.x =element_text(size=15,face="bold", color = "black"),
            axis.text.y =element_text(size=15,face="bold", color = "black"),
            axis.title=element_text(size=15,face="bold", color = "black"),
            aspect.ratio = 1
      )
    ggsave(paste0("/Users/saeko/Unmeasured/plot/comp_year/comp_SRX_year_", tgt_year, ".pdf"), p1, width = 9, height = 7)
  }
  
}