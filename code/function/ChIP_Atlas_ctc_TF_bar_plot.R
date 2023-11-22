ChIP_Atlas_ctc_TF_bar_plot <- function(calc_opt){
  
  # version1: FPKM全TFの上位25% (2.77)
  library(tidyverse)
  
  if(calc_opt == TRUE){
    
    # ChIP-Atlas --------
    df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/SRX_date_20231004.tsv")
    df2 <- df %>% filter(Antigen != "Epitope tags" & Antigen != "GFP" & Antigen != "Biotin") %>% 
      group_by(Cell_type_class, Antigen) %>% summarize(num = n())
    df3 <- df2 %>% pivot_wider(names_from = Cell_type_class, values_from = num)
    row_lab <- df3$Antigen
    df4 <- df3 %>% select(-Antigen)
    df5 <- df4 %>% as.matrix()
    rownames(df5) <- row_lab
    df5[which(is.na(df5))] <- 0
    saveRDS(df5, "/Users/saeko/Unmeasured/data/TF_CTC_sample_num_new_20231004.rds")
    
    # ChIP-Atlas + Refseq df ----
    mt_Chip <- readRDS("/Users/saeko/Unmeasured/data/TF_CTC_sample_num_new_20231004.rds")
    mt_RNA <- readRDS("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/RNAseq/RefEx_expression_RNAseq10_human_PRJEB2445_heat.rds")
    share_TF <- intersect(rownames(mt_RNA), rownames(mt_Chip)) 
    
    mt_RNA <- mt_RNA[share_TF, ]
    mt_Chip <- mt_Chip[share_TF, ]
    
    # Antigenの行を2つの行列で揃える
    order_RNA<- order(rownames(mt_RNA))
    order_Chip<- order(rownames(mt_Chip))
    mt_RNA2 <- mt_RNA[order_RNA,]
    mt_Chip2 <-mt_Chip[order_Chip,]
    
    # Ctc, tissueを2つの行列で揃える-----
    colnames(mt_Chip2)
    colnames(mt_RNA2)
    mt_Chip3 <- mt_Chip2 %>% as_tibble() %>% select(-Bone, -Embryo, -Epidermis, -`No description`, -Others, -Pancreas, -Placenta, -`Pluripotent stem cell`, -Unclassified)
    
    # ChIP-AtlasのCell type classアノテーションを、RefEXに合わせる -----
    mt_Chip4 <- mt_Chip3 %>% mutate(muscular = rowSums(mt_Chip3[, c("Cardiovascular", "Muscle")])) %>% select(-Cardiovascular, -Muscle) #CardiovascularとMuscleはmuscular"へ
    mt_Chip5 <- mt_Chip4 %>% mutate(reproductive = rowSums(mt_Chip3[, c("Gonad", "Prostate", "Uterus")])) %>% select(-Gonad, -Prostate, -Uterus) #GonadとProstateとUterusは"reproductive"へ
    
    mt_Chip6 <- mt_Chip5 %>% select("Neural",  "Blood", "Adipocyte", "reproductive", "muscular", "Digestive tract", "Liver", "Lung", "Kidney", "Breast") 
    colnames(mt_Chip6) <- colnames(mt_RNA2)
    mt_Chip6 <- mt_Chip6 %>% as.matrix()
    rownames(mt_Chip6) <- rownames(mt_Chip2)
    saveRDS(mt_RNA2, "/Users/saeko/Unmeasured/data/mt_RNA2.rds")
    saveRDS(mt_Chip6, "/Users/saeko/Unmeasured/data/mt_Chip6.rds")
    
  }
  
  mt_RNA2 <- readRDS("/Users/saeko/Unmeasured/data/mt_RNA2.rds")
  mt_Chip6 <- readRDS("/Users/saeko/Unmeasured/data/mt_Chip6.rds")
  
  # FPKM上位25%かつChIP-seqサンプルが0であるTF-CTCの組み合わせを1(=Unmeasured), それ以外を0にしたテーブルを作成
  mt_ovl <- matrix(nrow = nrow(mt_Chip6), ncol = ncol(mt_Chip6))
  mt_RNA3 <- matrix(nrow = nrow(mt_Chip6), ncol = ncol(mt_Chip6))
  for (i in 1:nrow(mt_Chip6)) {
    for (j in 1:ncol(mt_Chip6)) {
      tgt_Chip <- mt_Chip6[i,j]
      tgt_RNA <- mt_RNA2[i,j]
      
      if(tgt_RNA >= 2.77){
        mt_RNA3[i,j] <- 1
      }else{
        mt_RNA3[i,j] <- 0
      }
      
      if(tgt_Chip == 0 & tgt_RNA >= 2.77){
        mt_ovl[i,j] <- 1
      }else{
        mt_ovl[i,j] <- 0
      }
    }
  }
  rownames(mt_RNA3) <- rownames(mt_Chip6)
  colnames(mt_RNA3) <- colnames(mt_Chip6)
  rownames(mt_ovl) <- rownames(mt_Chip6)
  colnames(mt_ovl) <- colnames(mt_Chip6)
  saveRDS(mt_RNA3, "/Users/saeko/Unmeasured/data/mt_RNA3_FPKM_thre277.rds") 
  saveRDS(mt_ovl, "/Users/saeko/Unmeasured/data/mt_ovl.rds")
  
  mt_RNA3 <- readRDS("/Users/saeko/Unmeasured/data/mt_RNA3_FPKM_thre277.rds")
  mt_ovl <- readRDS("/Users/saeko/Unmeasured/data/mt_ovl.rds")
  
  column_sums_RNA3 <- colSums(mt_RNA3[, 1:ncol(mt_RNA3)]) #FPKMで定義されたTF-CTCの数
  column_sums_ovl <- colSums(mt_ovl[, 1:ncol(mt_ovl)]) #FPKMかつChIP-seqの有無で定義したUnmeasuredのTF-CTCの数
  
  df_plot1 <- tibble(Cell_type_class = names(column_sums_RNA3), 
                     TF_CTC_FPKM_num = column_sums_RNA3,
                     Unmeasured = column_sums_ovl)
  df_plot2 <- df_plot1 %>% mutate(Unmeasured_percentage = (Unmeasured/TF_CTC_FPKM_num)*100) %>%
    mutate(Measured = TF_CTC_FPKM_num - Unmeasured) %>% 
    pivot_longer(-c(Cell_type_class, TF_CTC_FPKM_num, Unmeasured_percentage), names_to = "measure", values_to = "number")
  
  p1 <- df_plot2 %>% 
    ggplot(aes(x = reorder(Cell_type_class, -number), y = number, fill = measure)) +
    geom_bar(stat = "identity", width = 0.7) +
    xlab("Cell type class (tissue)")+
    ylab("Number of expressed TFs")+
    scale_fill_manual(values = c("gray", "blue3")) +
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
          aspect.ratio = 0.5
    )
  
  # unmeasured ratio 
  df_plot3 <- df_plot1 %>% mutate(Unmeasured_percentage = (Unmeasured/TF_CTC_FPKM_num)*100)
  p2 <- df_plot3 %>% 
    ggplot(aes(x = reorder(Cell_type_class, -Unmeasured_percentage), y = Unmeasured_percentage)) +
    geom_point(size = 3) +
    ylim(c(0, 100)) +
    xlab("Cell type class (tissue)")+
    ylab("Unmeasured (%)")+
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
          aspect.ratio = 0.5
    )
  
  return(list(p1, p2))
}