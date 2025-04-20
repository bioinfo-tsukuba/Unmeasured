GWAS_cover_ratio_year_plot <- function(){
  
  library(tidyverse)
  library(ggrepel)
  library(RColorBrewer)

  ctc_plot_df2 <- read_tsv("/Users/saeko/Unmeasured/data/GWAS_cover_ratio/ctc_plot_df2.tsv")
  color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG"))
  ctc_plot_df3 <- ctc_plot_df2 %>% mutate(label = ifelse(year == 2023, Cell_type_class, ""))
  
  p1 <- ctc_plot_df3 %>% drop_na(Cell_type_class) %>% 
    filter(Cell_type_class != "Unclassified" & Cell_type_class != "Others") %>%
    ggplot(aes(x = year, y = GWAS_cover_ratio, color = Cell_type_class, label = label))  +
    geom_point() +
    geom_line()+
    geom_label_repel(size = 10) +
    scale_color_manual(values = color_list) +
    labs(
      x = "Year",
      y = "GWAS-SNP cover ratio",
      colour = "Cell type class"
    )+
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
  
  p2 <- ctc_plot_df3 %>% drop_na(Cell_type_class) %>%
    filter(Cell_type_class != "Unclassified" & Cell_type_class != "Others") %>%
    ggplot(aes(x = cum_uniq_TF, y = GWAS_cover_ratio, color = Cell_type_class, label = label))  +
    geom_point() +
    geom_line()+
    geom_label_repel(size = 10) +
    scale_color_manual(values = color_list) +
    labs(
      x = "Number of unique TFs",
      y = "GWAS-SNP cover ratio",
      colour = "Gears"
    )+
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
  
  
  return(list(p1, p2))
  
}