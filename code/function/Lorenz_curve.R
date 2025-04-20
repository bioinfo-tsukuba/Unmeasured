Lorenz_curve <- function(){
  
  library(ineq) 
  library(tidyverse)
  
  remap2 <- read_tsv("/Users/saeko/Unmeasured/data/ReMap2022/Remap_meta_human2.tsv")
  GTRD_meta_human2 <- read_tsv("/Users/saeko/Unmeasured/data/GTRD/GTRD_meta_human2.tsv")
  ChIP_annotation <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/SRX_date_20231004.tsv")
  
  # per TF -----
  ## remap2
  remap2_TF_num <- remap2 %>% group_by(TF) %>% summarise(n = n_distinct(number))
  remap2_TF_num2 <- remap2_TF_num %>%
    arrange(n) %>%  # 値が小さい順にソート
    mutate(
      cumulative_ChIPseq = cumsum(n) / sum(n),  # 累積シェア（y軸）
      cumulative_TF = (row_number()) / n(),  # 累積TF割合（x軸）
      db = "remap2"
    )
  
  p_TF_remap2 <- remap2_TF_num2 %>%
    ggplot(aes(x = cumulative_TF, y = cumulative_ChIPseq)) +
    geom_line(color = "black", size = 1.2) +  # ローレンツ曲線
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # 完全平等線
    scale_linetype_manual(values = c("Lorenz Curve" = "solid", "Equality Line" = "dashed")) +  # 線の種類
    labs(title = "Lorenz curve (remap2)",
         x = "Cumulative normalized rank of TF", 
         y = "Comulative normalized ChIP-seq")+
    annotate("text", x = 0.2, y = 0.65, label = "Lorenz_curve", color = "black", size = 5, hjust = 0) +
    annotate("text", x = 0.2, y = 0.7, label = "Line of perfect equality", color = "black", size = 5, hjust = 0) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.65, yend = 0.65, color = "black", linetype = "solid", size = 1) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.7, yend = 0.7, color = "black", linetype = "dashed", size = 1)+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=15, color = "black"),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=15, color = "black"),
          axis.title=element_text(size=15, color = "black"),
          aspect.ratio = 1
    )
    
  
  ## GTRD
  GTRD_TF_num <- GTRD_meta_human2 %>% group_by(gene_symbol) %>% summarise(n = n_distinct(id))
  GTRD_TF_num2 <- GTRD_TF_num %>%
    arrange(n) %>%  # 値が小さい順にソート
    mutate(
      cumulative_ChIPseq = cumsum(n) / sum(n),  # 累積シェア（y軸）
      cumulative_TF = (row_number()) / n(),  # 累積TF割合（x軸）
      db = "GTRD"
    )
  
  p_TF_GTRD <- GTRD_TF_num2 %>%
    ggplot(aes(x = cumulative_TF, y = cumulative_ChIPseq)) +
    geom_line(color = "black", size = 1.2) +  # ローレンツ曲線
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # 完全平等線
    scale_linetype_manual(values = c("Lorenz Curve" = "solid", "Equality Line" = "dashed")) +  # 線の種類
    labs(title = "Lorenz curve (GTRD)",
         x = "Cumulative normalized rank of TF", 
         y = "Comulative normalized ChIP-seq")+
    annotate("text", x = 0.2, y = 0.65, label = "Lorenz_curve", color = "black", size = 5, hjust = 0) +
    annotate("text", x = 0.2, y = 0.7, label = "Line of perfect equality", color = "black", size = 5, hjust = 0) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.65, yend = 0.65, color = "black", linetype = "solid", size = 1) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.7, yend = 0.7, color = "black", linetype = "dashed", size = 1)+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=15, color = "black"),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=15, color = "black"),
          axis.title=element_text(size=15, color = "black"),
          aspect.ratio = 1
    )
  
  ## ChIP-Atlas
  ChIP_TF_num <- ChIP_annotation %>% group_by(Antigen) %>% summarise(n = n_distinct(SRX))
  ChIP_TF_num2 <- ChIP_TF_num %>%
    arrange(n) %>%  # 値が小さい順にソート
    mutate(
      cumulative_ChIPseq = cumsum(n) / sum(n),  # 累積シェア（y軸）
      cumulative_TF = (row_number()) / n(),  # 累積TF割合（x軸）
      db = "ChIPatlas"
    ) 
  
  p_TF_ChIP <- ChIP_TF_num2 %>%
    ggplot(aes(x = cumulative_TF, y = cumulative_ChIPseq)) +
    geom_line(color = "black", size = 1.2) +  # ローレンツ曲線
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # 完全平等線
    scale_linetype_manual(values = c("Lorenz Curve" = "solid", "Equality Line" = "dashed")) +  # 線の種類
    labs(title = "Lorenz curve (ChIP-ATtlas)",
         x = "Cumulative normalized rank of TF", 
         y = "Comulative normalized ChIP-seq")+
    annotate("text", x = 0.2, y = 0.65, label = "Lorenz_curve", color = "black", size = 5, hjust = 0) +
    annotate("text", x = 0.2, y = 0.7, label = "Line of perfect equality", color = "black", size = 5, hjust = 0) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.65, yend = 0.65, color = "black", linetype = "solid", size = 1) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.7, yend = 0.7, color = "black", linetype = "dashed", size = 1)+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=15, color = "black"),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=15, color = "black"),
          axis.title=element_text(size=15, color = "black"),
          aspect.ratio = 1
    )
  
  ## combined graph
  colnames(GTRD_TF_num2) <- colnames(remap2_TF_num2)
  colnames(ChIP_TF_num2) <- colnames(remap2_TF_num2)
  df_TF_all <- remap2_TF_num2 %>% add_row(GTRD_TF_num2) %>% add_row(ChIP_TF_num2)
  
  p_TF_all <- df_TF_all %>%
    ggplot(aes(x = cumulative_TF, y = cumulative_ChIPseq, color = db)) +
    geom_line( size = 1.2) +  # ローレンツ曲線
    scale_color_manual(values =  c("red", "blue", "green4")) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # 完全平等線
    scale_linetype_manual(values = c("Lorenz Curve" = "solid", "Equality Line" = "dashed")) +  # 線の種類
    labs(title = "Lorenz curve",
         x = "Cumulative normalized rank of TF", 
         y = "Comulative normalized ChIP-seq")+
    annotate("text", x = 0.2, y = 0.6, label = "ChIP-Atlas", color = "red", size = 5, hjust = 0) +
    annotate("text", x = 0.2, y = 0.7, label = "GTRD", color = "blue", size = 5, hjust = 0) +
    annotate("text", x = 0.2, y = 0.8, label = "remap2", color = "green4", size = 5, hjust = 0) +
    annotate("text", x = 0.2, y = 0.9, label = "Line of perfect equality", color = "black", size = 5, hjust = 0) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.6, yend = 0.6, color = "red", linetype = "solid", size = 1) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.7, yend = 0.7, color = "blue", linetype = "solid", size = 1) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.8, yend = 0.8, color = "green4", linetype = "solid", size = 1) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.9, yend = 0.9, color = "black", linetype = "dashed", size = 1)+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=15, color = "black"),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=15, color = "black"),
          axis.title=element_text(size=15, color = "black"),
          aspect.ratio = 1
    )
  
  # per Cell type -----
  ## remap2
  remap2_Ct_num <- remap2 %>% group_by(Cell_type) %>% summarise(n = n_distinct(number))
  remap2_Ct_num2 <- remap2_Ct_num %>%
    arrange(n) %>%  # 値が小さい順にソート
    mutate(
      cumulative_ChIPseq = cumsum(n) / sum(n),  # 累積シェア（y軸）
      cumulative_Ct = (row_number()) / n(),  # 累積Cell type割合（x軸）
      db = "remap2"
    )
  
  p_Ct_remap2 <- remap2_Ct_num2 %>%
    ggplot(aes(x = cumulative_Ct, y = cumulative_ChIPseq)) +
    geom_line(color = "black", size = 1.2) +  # ローレンツ曲線
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # 完全平等線
    scale_linetype_manual(values = c("Lorenz Curve" = "solid", "Equality Line" = "dashed")) +  # 線の種類
    labs(title = "Lorenz curve (remap2)",
         x = "Cumulative normalized rank of Cell type", 
         y = "Comulative normalized ChIP-seq")+
    annotate("text", x = 0.2, y = 0.65, label = "Lorenz_curve", color = "black", size = 5, hjust = 0) +
    annotate("text", x = 0.2, y = 0.7, label = "Line of perfect equality", color = "black", size = 5, hjust = 0) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.65, yend = 0.65, color = "black", linetype = "solid", size = 1) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.7, yend = 0.7, color = "black", linetype = "dashed", size = 1)+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=15, color = "black"),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=15, color = "black"),
          axis.title=element_text(size=15, color = "black"),
          aspect.ratio = 1
    )
  
  
  ## GTRD
  GTRD_Ct_num <- GTRD_meta_human2 %>% group_by(title) %>% summarise(n = n_distinct(id))
  GTRD_Ct_num2 <- GTRD_Ct_num %>%
    arrange(n) %>%  # 値が小さい順にソート
    mutate(
      cumulative_ChIPseq = cumsum(n) / sum(n),  # 累積シェア（y軸）
      cumulative_Ct = (row_number()) / n(),  # 累積Cell type割合（x軸）
      db = "GTRD"
    )
  
  p_Ct_GTRD <- GTRD_Ct_num2 %>%
    ggplot(aes(x = cumulative_Ct, y = cumulative_ChIPseq)) +
    geom_line(color = "black", size = 1.2) +  # ローレンツ曲線
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # 完全平等線
    scale_linetype_manual(values = c("Lorenz Curve" = "solid", "Equality Line" = "dashed")) +  # 線の種類
    labs(title = "Lorenz curve (GTRD)",
         x = "Cumulative normalized rank of Cell type", 
         y = "Comulative normalized ChIP-seq")+
    annotate("text", x = 0.2, y = 0.65, label = "Lorenz_curve", color = "black", size = 5, hjust = 0) +
    annotate("text", x = 0.2, y = 0.7, label = "Line of perfect equality", color = "black", size = 5, hjust = 0) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.65, yend = 0.65, color = "black", linetype = "solid", size = 1) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.7, yend = 0.7, color = "black", linetype = "dashed", size = 1)+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=15, color = "black"),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=15, color = "black"),
          axis.title=element_text(size=15, color = "black"),
          aspect.ratio = 1
    )
  
  
  
  ## ChIP-Atlas
  ChIP_Ct_num <- ChIP_annotation %>% group_by(Cell_type) %>% summarise(n = n_distinct(SRX))
  ChIP_Ct_num2 <- ChIP_Ct_num %>%
    arrange(n) %>%  # 値が小さい順にソート
    mutate(
      cumulative_ChIPseq = cumsum(n) / sum(n),  # 累積シェア（y軸）
      cumulative_Ct = (row_number()) / n(),  # 累積Cell type割合（x軸）
      db = "ChIPatlas"
    )
  
  p_Ct_ChIP <- ChIP_Ct_num2 %>%
    ggplot(aes(x = cumulative_Ct, y = cumulative_ChIPseq)) +
    geom_line(color = "black", size = 1.2) +  # ローレンツ曲線
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # 完全平等線
    scale_linetype_manual(values = c("Lorenz Curve" = "solid", "Equality Line" = "dashed")) +  # 線の種類
    labs(title = "Lorenz curve (ChIP-Atlas)",
         x = "Cumulative normalized rank of Cell type", 
         y = "Comulative normalized ChIP-seq")+
    annotate("text", x = 0.2, y = 0.65, label = "Lorenz_curve", color = "black", size = 5, hjust = 0) +
    annotate("text", x = 0.2, y = 0.7, label = "Line of perfect equality", color = "black", size = 5, hjust = 0) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.65, yend = 0.65, color = "black", linetype = "solid", size = 1) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.7, yend = 0.7, color = "black", linetype = "dashed", size = 1)+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=15, color = "black"),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=15, color = "black"),
          axis.title=element_text(size=15, color = "black"),
          aspect.ratio = 1
    )
  
  ## combined graph
  colnames(GTRD_Ct_num2) <- colnames(remap2_Ct_num2)
  colnames(ChIP_Ct_num2) <- colnames(remap2_Ct_num2)
  df_Ct_all <- remap2_Ct_num2 %>% add_row(GTRD_Ct_num2) %>% add_row(ChIP_Ct_num2)
  
  p_Ct_all <- df_Ct_all %>%
    ggplot(aes(x = cumulative_Ct, y = cumulative_ChIPseq, color = db)) +
    geom_line( size = 1.2) +  # ローレンツ曲線
    scale_color_manual(values =  c("red", "blue", "green4")) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # 完全平等線
    scale_linetype_manual(values = c("Lorenz Curve" = "solid", "Equality Line" = "dashed")) +  # 線の種類
    labs(title = "Lorenz curve",
         x = "Cumulative normalized rank of Cell type", 
         y = "Comulative normalized ChIP-seq")+
    annotate("text", x = 0.2, y = 0.6, label = "ChIP-Atlas", color = "red", size = 5, hjust = 0) +
    annotate("text", x = 0.2, y = 0.7, label = "GTRD", color = "blue", size = 5, hjust = 0) +
    annotate("text", x = 0.2, y = 0.8, label = "remap2", color = "green4", size = 5, hjust = 0) +
    annotate("text", x = 0.2, y = 0.9, label = "Line of perfect equality", color = "black", size = 5, hjust = 0) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.6, yend = 0.6, color = "red", linetype = "solid", size = 1) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.7, yend = 0.7, color = "blue", linetype = "solid", size = 1) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.8, yend = 0.8, color = "green4", linetype = "solid", size = 1) +
    annotate("segment", x = 0.05, xend = 0.18, y = 0.9, yend = 0.9, color = "black", linetype = "dashed", size = 1)+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=15, color = "black"),
          axis.text.x =element_blank(),
          axis.text.y =element_text(size=15, color = "black"),
          axis.title=element_text(size=15, color = "black"),
          aspect.ratio = 1
    )
  
  
  return(list(p_TF_remap2, p_TF_GTRD, p_TF_ChIP, p_TF_all, p_Ct_remap2, p_Ct_GTRD, p_Ct_ChIP, p_Ct_all))
  
}