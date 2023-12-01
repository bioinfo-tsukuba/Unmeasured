combine_year <- function(tgt_ct, tgt_ct2, tgt_region, tgt_threshold){
  
  library(tidyverse)
  library(RColorBrewer)
  library(patchwork)
  
  print(tgt_ct)
  print(tgt_ct2)
  
  RR_cover_ratio_tib <- suppressMessages(read_tsv(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/RR_ratio_table/RR_cover_ratio_", tgt_ct2, "_",tgt_region, ".tsv")))
  
  # expressed geneのうち、thresholdの数以上regulatory TFがあるgeneの割合
  #color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG"))
  color_list <- c(brewer.pal(10,"Spectral"))
  p1 <- RR_cover_ratio_tib %>% 
    filter(threshold == tgt_threshold) %>%
    ggplot(aes(x = year, y = RR_cover_ratio)) +
    geom_point() +
    geom_line() +
    ylim(c(0,1))+
    xlab("Year")+
    ylab("Reg-TF cover ratio")+
    labs(color = "Threshold")+
    #scale_color_manual(values = color_list) +
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
          axis.text=element_text(size=15,face="bold", color = "black"),
          axis.text.x =element_text(size=15,face="bold", color = "black", angle = 45, hjust = 1),
          axis.text.y =element_text(size=15,face="bold", color = "black"),
          axis.title=element_text(size=15,face="bold", color = "black"),
          aspect.ratio = 0.5
    )
  
  uniq_TF_df <- suppressMessages(read_tsv(paste0("/Users/saeko/Unmeasured/data/uniqTF_table/", tgt_ct2, "_", tgt_region, "bp.tsv")))
  p2 <- uniq_TF_df %>% ggplot(aes(x = year, y = uniq_TF_num)) +
    geom_point()+
    geom_line()+
    xlab("Year")+
    ylab("Number of unique TFs")+
    ggtitle(tgt_ct2)+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = "none",
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
  
  SRX_df <- suppressMessages(read_tsv(paste0("/Users/saeko/Unmeasured/data/SRX_table/", tgt_ct2, "_", tgt_region,  "bp.tsv")))
  p3 <- SRX_df %>% ggplot(aes(x = year, y = SRX_num)) +
    geom_point()+
    geom_line()+
    xlab("Year")+
    ylab("Number of ChIP-seq")+
    ggtitle(tgt_ct2)+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = "none",
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
  
  system(paste0("mkdir -p /Users/saeko/Unmeasured/plot/combine_year/threshold_", tgt_threshold))
  pdf(file = paste0("/Users/saeko/Unmeasured/plot/combine_year/threshold_", tgt_threshold, "/", tgt_ct2,  "_threshold_", tgt_threshold,"_", tgt_region, "bp.pdf" ), width = 5, height = 9)
  print(p1 / p2 / p3)
  dev.off()
  
}
