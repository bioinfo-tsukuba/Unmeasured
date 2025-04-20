Dotplot_fam_ctc_unmeasured_v2 <- function(calc_opt, tgt_threshold){
  
  # 「１度もChIP-seqされていないTF」も含めたバージョン
  
  library(tidyverse)
  library(ggside)
  
  # ChIP-Atlas
  df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/SRX_date_20231004.tsv")
  df <- df %>% arrange(year, month, day, time) %>% mutate(order = 1:nrow(df)) %>% 
    filter(Cell_type_class != "Unclassified" & Cell_type_class != "No description") %>%
    filter(Antigen != "GFP" & Antigen != "Epitope tags")
  total_sample <- length(unique(df$SRX))
  ChIPatlas_TF <- df$Antigen %>% unique()
  
  df2 <- df %>% select(-c(year, month, day)) %>%
    group_by(Antigen, Cell_type_class) %>%
    summarise(n_sample = n())
  
  # unmeasured TF
  df_S1 <- read_tsv("/Users/saeko/Unmeasured/data/Human_transcriptional_factors_Cell_2018/Supple_table_S1.tsv")
  df_S1_TF <- df_S1 %>% filter(IsTF == "Yes")
  length(unique(df_S1_TF$TF))
  unmeasured_TF <- setdiff(unique(df_S1_TF$TF), ChIPatlas_TF)
  length(unmeasured_TF)
  
  df2_join <- tibble(Antigen = unmeasured_TF, Cell_type_class = "Blood", n_sample = 0)
  df2 <- rbind(df2, df2_join)
  
  df_fam <- readRDS("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS_paper_public/data/Fig2/obj/df_fam.rds")
  tgt_IDs <- df$SRX %>% unique() %>% as.character()
  tgt_df_fam <- df_fam %>% filter(ID %in% tgt_IDs) %>% select(Antigen, Family) %>% distinct()
  
  df3 <- df2 %>% left_join(tgt_df_fam, by = "Antigen") %>%
    filter(Family != "Unknown" & Family != "NA" & Family != "No_annotation")
  df3_join <- df3 %>% select(Antigen, Family) %>% distinct()
  
  # Unmeasured(%)の計算
  if(calc_opt == TRUE){
    
    # RNAの発現に応じてtotal TF tableを作成
    mt_RNA4 <- readRDS("/Users/saeko/Unmeasured/data/mt_RNA4.rds")
    mt_Chip7 <- readRDS("/Users/saeko/Unmeasured/data/mt_Chip7.rds")
    
    # FPKM細胞ごとに上位25% or 50%であるTF-CTCの組み合わせをcount
    mt_totalTF <- matrix(nrow = nrow(mt_Chip6), ncol = ncol(mt_Chip6))
    for (i in 1:nrow(mt_Chip6)) {
      for (j in 1:ncol(mt_Chip6)) {
        
        tgt_TF <- rownames(mt_Chip6)[i]
        tgt_tissue <- colnames(mt_Chip6)[j]
        threshold <- quantile(mt_RNA4[, tgt_tissue], tgt_threshold) 
        #threshold <- quantile(mt_RNA2[, j], tgt_threshold) #第三四分位数を求める
        if(tgt_TF %in% rownames(mt_RNA4)){
          tgt_RNA <- mt_RNA4[tgt_TF,tgt_tissue]
          
          if(tgt_RNA >= threshold){
            mt_totalTF[i,j] <- 1
          }else{
            mt_totalTF[i,j] <- 0
          }
        }
      }
    }
    rownames(mt_totalTF) <- rownames(mt_Chip7)
    colnames(mt_totalTF) <- colnames(mt_Chip7)
    saveRDS(mt_totalTF, paste0("/Users/saeko/Unmeasured/data/mt_total_expressed_TF_",tgt_threshold,"_v2.rds"))
  }
  
  mt_totalTF <- readRDS(paste0("/Users/saeko/Unmeasured/data/mt_total_expressed_TF_", tgt_threshold,"_v2.rds"))
  rownames_TF <- rownames(mt_totalTF)
  mt_totalTF_tib <- mt_totalTF %>% as_tibble() %>% mutate(Antigen = rownames_TF) %>%
    left_join(df3_join, by = "Antigen") %>%
    pivot_longer(-c(Antigen, Family), names_to = "Cell_type_class", values_to = "expressed_TF") %>%
    filter(is.na(Family) == FALSE & Family != "No_annotation") %>%
    group_by(Family, Cell_type_class) %>%
    summarise(total_expressed_TF = sum(expressed_TF))
  
  mt_Chip7 <- readRDS("/Users/saeko/Unmeasured/data/mt_Chip7.rds")
  rownames_TF_chip <- rownames(mt_Chip6)
  mt_Chip7_tib <- mt_Chip7 %>% as_tibble() %>% mutate(Antigen = rownames_TF_chip) %>%
    left_join(df3_join, by = "Antigen") %>%
    pivot_longer(-c(Antigen, Family), names_to = "Cell_type_class", values_to = "num_ChIPseq") %>%
    filter(is.na(Family) == FALSE & Family != "No_annotation") %>%
    mutate(uniq_measured_TF = ifelse(num_ChIPseq >= 1, 1, 0)) %>%
    group_by(Family, Cell_type_class) %>%
    summarise(total_measured_uniq_TF = sum(uniq_measured_TF))
  
  tib_plot <-  mt_totalTF_tib %>% left_join(mt_Chip7_tib, by = c("Family", "Cell_type_class")) %>%
    mutate(ratio_unmeasured = ifelse(total_measured_uniq_TF < total_expressed_TF & total_expressed_TF != 0, 
                                     1-total_measured_uniq_TF/total_expressed_TF, 0),
           unmeasured_TF_num = total_expressed_TF - total_measured_uniq_TF) %>%
    mutate(unmeasured_TF_num = ifelse(unmeasured_TF_num<0, 0, unmeasured_TF_num))
  
  tmp <- tib_plot
  tmp2 <- tmp %>% group_by(Cell_type_class) %>% summarise(n_sum_unmeasured_ctc = sum(unmeasured_TF_num))
  tmp3 <- tmp %>% left_join(tmp2, by = "Cell_type_class")
  
  tmp4 <- tmp %>% group_by(Family) %>% summarise(n_sum_unmeasured_fam = sum(unmeasured_TF_num))
  tmp5 <- tmp3 %>% left_join(tmp4, by = "Family")
  
  p <- tmp5 %>% 
    ggplot(aes(x = reorder(Family, -unmeasured_TF_num), y= reorder(Cell_type_class, unmeasured_TF_num), size = ratio_unmeasured, label = unmeasured_TF_num)) + 
    geom_point(aes(color = ratio_unmeasured)) +
    geom_text(color = "white", size = 2) +
    #scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,max(tib_plot$ratio_unmeasured)), oob = scales::squish) +
    scale_color_gradientn(
      #colours = c("blue", "white",  "red"),  # お好みの色ベクトル
      colours = c("white",  "red"),
      limits = c(0, max(tib_plot$ratio_unmeasured)), 
      oob = scales::squish
    )+
    xlab("TF family") +
    ylab("Cell type class") +
    theme(
      aspect.ratio = 0.8,
      plot.background = element_blank(),
      plot.title = element_text(size = unit(12, "pt")),
      plot.subtitle = element_text(size = unit(10, "pt")),
      plot.caption = element_text(size = unit(8, "pt")),
      plot.tag = element_text(size = unit(18, "pt")),
      axis.ticks = element_line(colour = "black", linewidth = unit(0.5, "pt")),
      axis.line = element_blank(), # element_line(linewidth = unit(0.5, "pt")),
      axis.title = element_text(size = unit(15, "pt")),
      axis.text.x = element_text(size= unit(15, "pt"), color = "black", angle = 90, hjust = 1),
      axis.text.y = element_text(size= unit(15, "pt"), color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = unit(0.25, "pt")),
      panel.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = unit(10, "pt")),
      legend.text = element_text(size = unit(10, "pt")),
      strip.text = element_text(colour = "black", size = unit(8, "pt")),
      strip.background = element_rect(fill = NA, color = NA))
  
  
  # p2 <- p + 
  #   geom_xsidecol(aes(x = Family, y = n_sample), data = tmp, stat = "identity", width = 0.5, fill = "blue4") +
  #   geom_ysidecol(aes(x = n_sample, y = Cell_type_class), data = tmp, stat = "identity", width = 0.5, fill = "blue4") +
  #   theme(
  #     legend.position = "none"
  #   ) 
  
  return(p)
  
}