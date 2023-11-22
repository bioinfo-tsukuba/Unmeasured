RR_cover_ratio_allTF <- function(tgt_ct, tgt_ct2, tgt_region){
  
  library(tidyverse)
  library(RColorBrewer)
  print(tgt_ct)
  print(tgt_ct2)
  
  RR_cover_ratio_tib <- read_tsv(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/RR_ratio_table/RR_cover_ratio_", tgt_ct2, "_",tgt_region, ".tsv"))
  
  # expressed geneのうち、thresholdの数以上regulatory TFがあるgeneの割合
  #color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG"))
  color_list <- c(brewer.pal(10,"Spectral"))
  p1 <- RR_cover_ratio_tib %>% 
    filter(threshold %in% 1:10) %>%
    ggplot(aes(x = year, y = RR_cover_ratio, color = factor(threshold))) +
    geom_point() +
    geom_line() +
    ylim(c(0,1))+
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
  
  return(list(p1))
}