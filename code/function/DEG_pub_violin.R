DEG_pub_violin <- function(){
  library(tidyverse)
  
  # /Users/saeko/Documents/MOCCS/important_chipseq_prediction/code/TFmarker_DEG_MeSH_v2_rev.Rより
  df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/DEG_TFmarker_pub/summary_allstats_v2.tsv")
  df <- df %>% drop_na(TF, Cell_type_class)
  df2 <- df %>% 
    mutate(label_marker = ifelse(is.na(Gene_type) == FALSE, "Marker", "No marker")) %>%
    mutate(label_chip_measure = ifelse(count_ChIPseq == 0, "Unmeasured", "Measured")) %>%
    select(-Gene_type) %>% distinct() 
  
  # Number of DEGs, violin plot -----
  p1 <- df2 %>% ggplot(aes(x = label_chip_measure, y = DEG_num_thre_p, fill = label_chip_measure)) +
    geom_violin() +
    geom_point() +
    scale_fill_manual(values = c("#D7191C", "blue4"))+
    xlab("Measure")+
    ylab("Number of DEGs")+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          #legend.position = "none",
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
  p2 <- df2 %>% ggplot(aes(x = label_chip_measure, y = log(pub_count, base = 10), fill = label_chip_measure)) +
    geom_violin() +
    geom_point() +
    scale_fill_manual(values = c("#D7191C", "blue4"))+
    xlab("Measure")+
    ylab("log10(Number of publication)")+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          #legend.position = "none",
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
  print("statistical test: DEGs")
  measured_DEG <- df2 %>% filter(label_chip_measure == "Measured") %>% .$DEG_num_thre_p
  unmeasured_DEG <- df2 %>% filter(label_chip_measure == "Unmeasured") %>% .$DEG_num_thre_p
  print(wilcox.test(measured_DEG, unmeasured_DEG, paired = FALSE))
  
  print("statistical test: publication")
  measured_pub <- df2 %>% filter(label_chip_measure == "Measured") %>% .$pub_count
  unmeasured_pub <- df2 %>% filter(label_chip_measure == "Unmeasured") %>% .$pub_count
  print(wilcox.test(measured_pub, unmeasured_pub, paired = FALSE))
  
  return(list(p1, p2))
}