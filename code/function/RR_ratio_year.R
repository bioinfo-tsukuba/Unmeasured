RR_cover_ratio <- function(tgt_ct, tgt_ct2, tgt_region){
  
  library(tidyverse)
  library(RColorBrewer)
  print(tgt_ct)
  print(tgt_ct2)
  
  # peak-overlap matrix -----
  df <- readRDS(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TSS_matrix/TSS_", tgt_region, "bp_overlaps_matrix2_pivot_loneger_annotated.rds")) # _selectedは、TF-Ctで1pairにpeak数maxで絞ったものなので、今回は使わない。
  df_wider <- df %>% select(ID, gene, peak_existing)%>% pivot_wider( names_from = gene, values_from = peak_existing)
  
  # geneをFPKM上位に絞る  ----
  if(file.exists(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/VirtualChIP/RNAseq/", tgt_ct, "_RNA.tsv"))==TRUE){
    RNA_df <- read_tsv(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/VirtualChIP/RNAseq/", tgt_ct, "_RNA.tsv"))
  }else{
    RNA_df <- read_tsv(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/VirtualChIP/RNAseq/CCLE/", tgt_ct, "_RNA.tsv"))
  }
  
  colnames(RNA_df) <- c("gene", "RPKM_or_TPM")
  
  gene <- RNA_df$gene
  RNA_df2 <- RNA_df %>% separate(col = gene, into = c("gene2", "gene3")) %>% mutate(gene = gene)
  tgt_gene <- colnames(df_wider)
  RNA_df3 <- RNA_df2 %>% filter(gene %in% tgt_gene | gene2 %in% tgt_gene)
  
  RPKMs <- RNA_df3 %>% arrange(desc(RPKM_or_TPM)) %>% .$RPKM_or_TPM %>% as.numeric()
  threshold <- RPKMs[0.25*nrow(RNA_df3)]
  expressed_genes <- RNA_df3 %>% filter(RPKM_or_TPM >= threshold & gene %in% tgt_gene) %>% .$gene
  print(paste0("Expressed genes in ", tgt_ct2, " : ", length(expressed_genes)))
  
  df_wider_selected <- df_wider[, c(expressed_genes, "ID")]
  dim(df_wider_selected)
  
  
  # yearごとのTF select df ------
  df_date <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/SRX_date_20231004.tsv")
  df_date <- df_date %>% arrange(year, month, day, time) %>% mutate(order = 1:nrow(df_date)) %>% 
    filter(Cell_type_class != "Unclassified" & Cell_type_class != "No description") %>%
    filter(Antigen != "GFP" & Antigen != "Epitope tags")
  df_date_selected <- df_date %>% filter(Cell_type == tgt_ct2)
  
  anno <- df_date_selected %>% as_tibble() %>% filter(SRX %in% df_wider_selected$ID) %>% distinct()
  colnames(anno) <- c("ID", colnames(anno)[2:11])
  df_wider_selected_annotated <- df_wider_selected %>% left_join(anno, by = "ID") %>%
    select(colnames(anno),  colnames(df_wider_selected)[1:(ncol(df_wider_selected)-1)]) %>%
    filter(ID %in% df_date_selected$SRX)
  write_tsv(df_wider_selected_annotated, paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/RR_ratio_table/df_wider_selected_annotated_", tgt_ct2, "_",tgt_region, "bp.tsv"))
  
  # make table for plot -----
  year_list <- unique(df_date$year) %>% sort()
  RR_cover_ratio_tib <- tibble()
  for (tgt_year in year_list) {
    print(tgt_year)
    for (tgt_threshold in 1:20) {
      
      # tgt yearまでのRR ratioを計算 (regTFがtgt_threshold以上あったらそのgeneはregTFありと判断)
      tgt_year_df <- df_wider_selected_annotated %>% filter(year %in% year_list[1]:tgt_year)
      tgt_year_df2 <- tgt_year_df %>% select(expressed_genes)
      
      # unique TFの数で絞りたい
      num_uniq_TF_list <- list()
      for (i in 1:length(expressed_genes)) {
        gene <- expressed_genes[i] %>% as.character()
        tgt_df <- tgt_year_df %>% 
          select("ID", "year", "month", "day", "time", "Antigen", "Cell_type_class", "Cell_type", gene) %>%
          filter(!!sym(gene) == 1)
        num_uniq_TF <- length(unique(tgt_df$Antigen)) %>% as.double()
        num_uniq_TF_list[[gene]] <- num_uniq_TF
      }
      colsums_df <- tibble(gene = names(num_uniq_TF_list), num_uniq_TF = as.double(num_uniq_TF_list)) %>% 
        mutate(year = tgt_year, Cell_type = tgt_ct2)
      #write_tsv(colsums_df, paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/RR_ratio_table/uniqTF_num_df_", tgt_ct2, "_",tgt_region, "bp.tsv"))
      regTF_exist_genes <- colsums_df %>% filter(num_uniq_TF > tgt_threshold) %>% .$gene
      tgt_RR_cover_ratio <- length(regTF_exist_genes) / length(expressed_genes)
      tgt_tib <- tibble(year = tgt_year, threshold = tgt_threshold, RR_cover_ratio = tgt_RR_cover_ratio)
      
      
      if(nrow(RR_cover_ratio_tib) == 0){
       RR_cover_ratio_tib <- tgt_tib
      }else{
        RR_cover_ratio_tib <- rbind(RR_cover_ratio_tib, tgt_tib)
      }
    }
  }
  RR_cover_ratio_tib$threshold <- RR_cover_ratio_tib$threshold %>% as.integer()
  write_tsv(RR_cover_ratio_tib, paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/RR_ratio_table/RR_cover_ratio_", tgt_ct2, "_",tgt_region, ".tsv"))
  
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
    ggtitle(paste0("Regulatory-TF cover ratio in expressed genes, ", tgt_ct2))+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          #legend.position = "none",
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
  
  return(list(p1))
}