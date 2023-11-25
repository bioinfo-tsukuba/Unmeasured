DEG_pub_violin <- function(){
  
  library(tidyverse)
  library(ggrepel)
  library(ggside)
  
  # /Users/saeko/Documents/MOCCS/important_chipseq_prediction/code/TFmarker_DEG_MeSH_v3_20231125.R より
  df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/DEG_TFmarker_pub/summary_allstats_v3.tsv")
  df <- df %>% drop_na(TF, Cell_type_class)
  df2 <- df %>% 
    mutate(label_marker = ifelse(is.na(Gene_type) == FALSE, "Marker", "No marker")) %>%
    mutate(label_chip_measure = ifelse(count_ChIPseq == 0, "Unmeasured", "Measured")) %>%
    select(-Gene_type) %>% distinct() 
  
  # Unmeasured vs Measured ----
  # Number of DEGs, violin plot -----
  p1 <- df2 %>% ggplot(aes(x = label_chip_measure, y = DEG_num_thre_p, fill = label_chip_measure)) +
    geom_violin() +
    geom_point() +
    scale_fill_manual(values = c("#D7191C", "blue4"))+
    xlab("Measure")+
    ylab("Number of DEGs")+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_line(color = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold", color = "black"),
          axis.text.x =element_text(size=10,face="bold", color = "black"),
          axis.text.y =element_text(size=10,face="bold", color = "black"),
          axis.title=element_text(size=15,face="bold", color = "black"),
          aspect.ratio = 2
    )
  
  # Number of publications, violin plot ----
  p2 <- df2 %>% 
    ggplot(aes(x = label_chip_measure, y = log(pub_count, base = 10), fill = label_chip_measure)) +
    geom_violin() +
    geom_point() +
    scale_fill_manual(values = c("#D7191C", "blue4"))+
    xlab("Measure")+
    ylab("log10(Number of publication)")+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_line(color = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold", color = "black"),
          axis.text.x =element_text(size=10,face="bold", color = "black"),
          axis.text.y =element_text(size=10,face="bold", color = "black"),
          axis.title=element_text(size=15,face="bold", color = "black"),
          aspect.ratio = 2
    )
  
  
  # 検定 ----
  print("statistical test: DEGs, Unmeasured vs Measured")
  measured_DEG <- df2 %>% filter(label_chip_measure == "Measured") %>% .$DEG_num_thre_p
  unmeasured_DEG <- df2 %>% filter(label_chip_measure == "Unmeasured") %>% .$DEG_num_thre_p
  print(wilcox.test(measured_DEG, unmeasured_DEG, paired = FALSE))
  
  print("statistical test: publication, Unmeasured vs Measured")
  measured_pub <- df2 %>% filter(label_chip_measure == "Measured") %>% .$pub_count
  unmeasured_pub <- df2 %>% filter(label_chip_measure == "Unmeasured") %>% .$pub_count
  print(wilcox.test(measured_pub, unmeasured_pub, paired = FALSE)) #paired = FALSEでwilcoxonの順位検定(not 符号付き)
  
  
  
  # TF marker vs No TF-marker ----
  # Number of DEGs, violin plot -----
  p3 <- df2 %>% ggplot(aes(x = label_marker, y = DEG_num_thre_p, fill = label_marker)) +
    geom_violin() +
    geom_point() +
    scale_fill_manual(values = c("#D7191C", "gray"))+
    xlab("TF marker")+
    ylab("Number of DEGs")+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_line(color = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold", color = "black"),
          axis.text.x =element_text(size=10,face="bold", color = "black"),
          axis.text.y =element_text(size=10,face="bold", color = "black"),
          axis.title=element_text(size=15,face="bold", color = "black"),
          aspect.ratio = 2
    )
  
  # Number of publications, violin plot ----
  p4 <- df2 %>% ggplot(aes(x = label_marker, y = log(pub_count, base = 10), fill = label_marker)) +
    geom_violin() +
    geom_point() +
    scale_fill_manual(values = c("#D7191C", "gray"))+
    xlab("TF marker")+
    ylab("log10(Number of publication)")+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_line(color = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold", color = "black"),
          axis.text.x =element_text(size=10,face="bold", color = "black"),
          axis.text.y =element_text(size=10,face="bold", color = "black"),
          axis.title=element_text(size=15,face="bold", color = "black"),
          aspect.ratio = 2
    )
  
  
  # 検定 ----
  print("statistical test: DEGs, TFmarker vs no TF-marker")
  marker_DEG <- df2 %>% filter(label_marker == "Marker") %>% .$DEG_num_thre_p
  no_marker_DEG <- df2 %>% filter(label_marker == "No marker") %>% .$DEG_num_thre_p
  print(wilcox.test(marker_DEG, no_marker_DEG, paired = FALSE))
  
  print("statistical test: publication, TFmarker vs no TF-marker")
  marker_pub <- df2 %>% filter(label_marker == "Marker") %>% .$pub_count
  no_marker_pub <- df2 %>% filter(label_marker == "No marker") %>% .$pub_count
  print(wilcox.test(marker_pub, no_marker_pub, paired = FALSE))
  
  
  # scatter plot -----
  p5 <- df2 %>% #select(-Gene_type) %>% 
    distinct() %>%
    ggplot(aes(x = log(pub_count, base = 10), y = DEG_num_thre_p, 
               color = label_chip_measure,
               label = key
    )) +
    geom_point() +
    geom_text_repel() +
    xlab("log10(Number of publication)")+
    ylab("Number of DEGs (corrected p value < 0.05)") +
    scale_color_manual(values = c("#D7191C", "blue4"))+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          #legend.position = "none",
          panel.grid.major = element_line(color = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold"),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=15,face="bold"),
          ggside.panel.scale = 0.2,
          aspect.ratio = 0.9
    )+
    geom_xsideviolin(aes(y = label_chip_measure), orientation = "y") + 
    geom_xsideboxplot(aes(y = label_chip_measure), orientation = "y") + 
    scale_xsidey_discrete(guide = guide_axis(angle = 0))+
    geom_ysideviolin(aes(x = label_chip_measure), orientation = "x") +
    geom_ysideboxplot(aes(x = label_chip_measure), orientation = "x") +
    scale_ysidex_discrete(guide = guide_axis(angle = 45))
  
  p6 <- df2 %>% #select(-Gene_type) %>% 
    distinct() %>%
    ggplot(aes(x = log(pub_count, base = 10), y = DEG_num_thre_p, 
               color = label_marker,
               label = key
    )) +
    geom_point() +
    geom_text_repel() +
    xlab("log10(Number of publication)")+
    ylab("Number of DEGs (corrected p value < 0.05)") +
    scale_color_manual(values = c("#D7191C", "gray"))+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          #legend.position = "none",
          panel.grid.major = element_line(color = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold"),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=15,face="bold"),
          ggside.panel.scale = 0.2,
          aspect.ratio = 0.9
    )+
    geom_xsideviolin(aes(y = label_marker), orientation = "y") + 
    geom_xsideboxplot(aes(y = label_marker), orientation = "y") + 
    scale_xsidey_discrete(guide = guide_axis(angle = 0))+
    geom_ysideviolin(aes(x = label_marker), orientation = "x") +
    geom_ysideboxplot(aes(x = label_marker), orientation = "x") +
    scale_ysidex_discrete(guide = guide_axis(angle = 45))
  
  
  # filter only TF marker ---
  p7 <- df2 %>% filter(label_marker == "Marker") %>%
    ggplot(aes(x = log(pub_count, base = 10), y = DEG_num_thre_p, 
               color = label_chip_measure,
               label = key
    )) +
    geom_point() +
    geom_text_repel() +
    ggtitle("TF markers")+ 
    xlab("log10(Number of publication)")+
    ylab("Number of DEGs (corrected p value < 0.05)") +
    scale_color_manual(values = c("gray", "blue4"))+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          #legend.position = "none",
          panel.grid.major = element_line(color = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold"),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=15,face="bold"),
          ggside.panel.scale = 0.2,
          aspect.ratio = 0.9
    )
  
  return(list(p1, p2, p3, p4, p5, p6, p7))
}