DEG_pub_violin_rev1 <- function(){
  
  # /Users/saeko/Unmeasured/code/analysis/Revise/ExpressedTF_def_ver2.Rmd がbase
  library(tidyverse)
  library(ggrepel)
  library(ggside)
  
  df_marker <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/DEG_TFmarker_pub/summary_allstats_v3.tsv")
  mt_unmeasured_tib <- readRDS("/Users/saeko/Unmeasured/data/tib_unmeasured_stepminer_Ctc.rds")
  unmeasured_join <- mt_unmeasured_tib %>% select(Antigen, Cell_type_class, num_unmeasured)
  colnames(unmeasured_join) <- c("TF", "Cell_type_class", "num_unmeasured")
  
  tib_marker_plot <- df_marker %>% 
    left_join(unmeasured_join, by = c("TF", "Cell_type_class")) %>%
    drop_na(TF, Cell_type_class) %>%
    mutate(label_marker = ifelse(is.na(Gene_type) == FALSE, "Marker", "No marker"),
           label_chip_measure = ifelse(num_unmeasured == 1, "Unmeasured", "Measured")) %>%
    filter(Cell_type_class != "Others" & is.na(label_chip_measure) == FALSE) %>%
    select(-Gene_type) %>% distinct() 
  
  p1 <- tib_marker_plot %>% ggplot(aes(x = label_chip_measure, y = DEG_num_thre_p, fill = label_chip_measure)) +
    geom_violin() +
    geom_point() +
    scale_fill_manual(values = c("gray", "blue4"))+
    xlab("Measure")+
    ylab("Number of DEGs")+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=30, color = "black"),
          axis.text.x =element_text(size=30, color = "black", angle = 45, hjust = 1),
          axis.text.y =element_text(size=30,color = "black"),
          axis.title.x =element_blank(), #element_text(size=30,face="bold", color = "black"),
          axis.title.y = element_text(size=30, color = "black"),
          aspect.ratio = 1
    )+
    stat_compare_means(
      method = "wilcox.test",         # t検定。Wilcoxonなら "wilcox.test"
      comparisons = list(c("Unmeasured", "Measured")),  # 比較したいグループのペア
      label = "p.signif",        # アスタリスク（*, **, ***）で表示
      size = 5                 # p値表示の文字サイズ
    )
  
  # Unmeasured vs Measured ----
  # Number of DEGs, violin plot -----
  
  p2 <- tib_marker_plot %>% 
    ggplot(aes(x = label_marker, y = DEG_num_thre_p, fill = label_marker)) +
    geom_violin() +
    geom_point() +
    scale_fill_manual(values = c("#D7191C", "gray"))+
    xlab("TF marker")+
    ylab("Number of DEGs")+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=30,color = "black"),
          axis.text.x =element_text(size=30, color = "black", angle = 45, hjust = 1),
          axis.text.y =element_text(size=30, color = "black"),
          axis.title.x =element_blank(), #element_text(size=30,face="bold", color = "black"),
          axis.title.y = element_text(size=30, color = "black"),
          aspect.ratio = 1
    )+
    stat_compare_means(
      method = "wilcox.test",         # t検定。Wilcoxonなら "wilcox.test"
      comparisons = list(c("Marker", "No marker")),  # 比較したいグループのペア
      label = "p.signif",        # アスタリスク（*, **, ***）で表示
      size = 3                 # p値表示の文字サイズ
    )
  
  p3 <- tib_marker_plot %>%
    distinct() %>%
    ggplot(aes(x = log(pub_count, base = 10), y = DEG_num_thre_p, 
               color = label_chip_measure,
               label = key
    )) +
    geom_point() +
    geom_text_repel() +
    ggtitle("TF-cell type class pairs in KnockTF") +
    xlab("log10(Number of publication)")+
    ylab("Number of DEGs") +
    scale_color_manual(values = c("#D7191C", "blue4"))+
    theme(plot.title = element_text(hjust = 0.5), 
          title = element_text(size=30,color = "black"), 
          legend.position = c(0.2,0.8),
          legend.title = element_text(size = 15, color = "black"),
          legend.text = element_text(size = 15, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=30, color = "black"),
          axis.text.x =element_text(size=30, color = "black"),
          axis.text.y =element_text(size=30,color = "black"),
          axis.title.x = element_text(size=30, color = "black"), 
          axis.title.y = element_text(size=30,color = "black"),
          aspect.ratio = 1
    )
  
  p4 <- tib_marker_plot %>%
    filter(label_marker == "Marker") %>%
    ggplot(aes(x = log(pub_count, base = 10), y = DEG_num_thre_p, 
               color = label_chip_measure,
               label = key
    )) +
    geom_point() +
    geom_text_repel() +
    ggtitle("TF markers")+ 
    xlab("log10(Number of publication)")+
    ylab("Number of DEGs") +
    scale_color_manual(values = c("gray", "blue4"))+
    geom_hline(yintercept = 1000, linetype = "dashed", color = "black")+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = c(0.2,0.9),
          legend.title = element_text(size = 15, color = "black"),
          legend.text = element_text(size = 15, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          title = element_text(size=30, color = "black"), 
          axis.text=element_text(size=30, color = "black"),
          axis.text.x =element_text(size=30, color = "black"),
          axis.text.y =element_text(size=30, color = "black"),
          axis.title.x = element_text(size=30, color = "black"),
          axis.title.y = element_text(size=30, color = "black"),
          aspect.ratio = 1
    )
  
  return(list(p1, p2, p3, p4))
}