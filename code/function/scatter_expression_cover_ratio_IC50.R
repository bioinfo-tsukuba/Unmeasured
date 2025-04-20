scatter_expression_cover_ratio_IC50 <- function(tgt_region){
  
  library(tidyverse)
  library(RColorBrewer)
  
  ct_list <- read_tsv("/Users/saeko/Unmeasured/data/cell_type_list_RNA_ChIPAtlas.txt", col_names = F)
  colnames(ct_list) <- c("ct", "ct2")
  
  for (tgt_threshold in 1:10) {
    
    print(tgt_threshold)
    df_binded1 <- c()
    df_binded_RNA <- c() #top25% expressed gene
    df_binded_RNA_row <- c()
    
    df <- readRDS(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TSS_matrix/TSS_", tgt_region, "bp_overlaps_matrix2_pivot_loneger_annotated.rds")) # _selectedは、TF-Ctで1pairにpeak数maxで絞ったものなので、今回は使わない。
    df_wider <- df %>% select(ID, gene, peak_existing)%>% pivot_wider( names_from = gene, values_from = peak_existing)
    
    for (tgt_row in 1:nrow(ct_list)) {
      tgt_ct <- ct_list[tgt_row, "ct"] %>% as.character() #RNA-seq
      tgt_ct2 <- ct_list[tgt_row, "ct2"] %>% as.character() #ChIP-Atlas
      
      # RNAseq
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
      
      #df_wider_selected <- df_wider[, c(expressed_genes, "ID")]
      #dim(df_wider_selected)
      
      RNA_df4 <- RNA_df3 %>% filter(gene %in% expressed_genes) %>% mutate(Cell_type = tgt_ct2)
      RNA_df3_2 <- RNA_df3 %>% mutate(Cell_type = tgt_ct2)
      
      if(tgt_row == 1){
        df_binded_RNA <- RNA_df4
        df_binded_RNA_row <- RNA_df3_2
      }else{
        df_binded_RNA <- rbind(df_binded_RNA, RNA_df4)
        df_binded_RNA_row <- rbind(df_binded_RNA_row, RNA_df3_2)
      }
      
      # UniqueTF, Ratio, scatter
      RR_cover_ratio_tib <- suppressMessages(read_tsv(paste0("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/RR_ratio_table/RR_cover_ratio_", tgt_ct2, "_",tgt_region, ".tsv")))
      uniq_TF_df <- suppressMessages(read_tsv(paste0("/Users/saeko/Unmeasured/data/uniqTF_table/", tgt_ct2, "_", tgt_region, "bp.tsv")))
      tgt_RR_cover_ratio_tib <- RR_cover_ratio_tib %>% filter(threshold == tgt_threshold)
      df_join1 <- tgt_RR_cover_ratio_tib %>% left_join(uniq_TF_df, by = "year") %>%
        mutate(Cell_type = tgt_ct2, tgt_region = tgt_region)
      
      if(tgt_row == 1){
        df_binded1 <- df_join1
      }else{
        df_binded1 <- rbind(df_binded1, df_join1)
      }
    }
    
    p1 <- df_binded_RNA  %>% ggplot(aes(x = reorder(Cell_type, -RPKM_or_TPM), y = log(RPKM_or_TPM), color = Cell_type)) +
      geom_violin() +
      geom_point(size = 0.5) +
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
    
    p2 <- df_binded_RNA  %>% filter(Cell_type %in% c("PANC-1", "PC-3", "HCT 116", "MG-63", "A2780", "HEK293-T-REx")) %>%
      ggplot(aes(x = reorder(Cell_type, -RPKM_or_TPM), y = log(RPKM_or_TPM), color = Cell_type)) +
      geom_violin() +
      geom_boxplot()+
      geom_point(size = 0.5) +
      ylab("log(TPM)") +
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
    
    p3 <- df_binded_RNA  %>% filter(!Cell_type %in% c("PANC-1", "PC-3", "HCT 116", "MG-63", "A2780", "HEK293-T-REx")) %>%
      ggplot(aes(x = reorder(Cell_type, -RPKM_or_TPM), y = log(RPKM_or_TPM), color = Cell_type)) +
      geom_violin() +
      geom_boxplot()+
      geom_point(size = 0.5) +
      ylab("log(RPKM)") +
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
    
    # 「遺伝子発現値の分布」を作れるデータから、発現遺伝子数（発現量 > 1 （例えば) である遺伝子数の数）を計算 
    df1 <- df_binded_RNA_row %>% 
      filter(Cell_type %in% c("PANC-1", "PC-3", "HCT 116", "MG-63", "A2780", "HEK293-T-REx")) %>%
      filter(RPKM_or_TPM > 3) %>% 
      group_by(Cell_type) %>% summarise(gene_num = n())
    df2 <- df_binded_RNA_row %>% 
      filter(!Cell_type %in% c("PANC-1", "PC-3", "HCT 116", "MG-63", "A2780", "HEK293-T-REx")) %>%
      filter(RPKM_or_TPM > 10) %>% 
      group_by(Cell_type) %>% summarise(gene_num = n())
    df3 <- rbind(df1, df2)
    
    # 「UniqueTF数が少なくてもRR ratioが上がっている度」をみれる図を作れるデータから、縦軸（Ratio (RR%) ）の値が 0.5 になるときの横軸（Unique TF ）の値 （IC50のようなイメージ）を cell typeごとに計算
    IC50_list <- list()
    for (tgt_ct2 in unique(df_binded1$Cell_type)) {
      tgt_df <- df_binded1 %>% filter(Cell_type == tgt_ct2) 
      x <- tgt_df$uniq_TF_num # x は 「Unique TF数とRRカバー率の散布図」の点の Unique TF (横軸）の値のベクトル
      y <- tgt_df$RR_cover_ratio # y は「Unique TF数とRRカバー率の散布図」の点の Ratio (RR%) （縦軸）の値のベクトル
      linear_interporation <- approx(x, y, n = 1000, ties = "ordered")
      target_y_value <- 0.5
      # linear_interporation には、線形補間をした結果として、縦軸の値が 0.5 になるときの横軸の値が格納される
      IC50_list[[tgt_ct2]] <-  linear_interporation$x[which.min(abs(linear_interporation$y - target_y_value))]
    }
    
    IC50_df <- tibble(Cell_type = names(IC50_list), IC50 = as.numeric(IC50_list))
    df4 <- IC50_df %>% left_join(df3, by = "Cell_type") %>% 
      mutate(Source = ifelse(Cell_type %in% c("PANC-1", "PC-3", "HCT 116", "MG-63", "A2780", "HEK293-T-REx"), "TPM", "RPKM")) 
    p4 <- df4 %>% ggplot(aes(x = gene_num, y = IC50, label = Cell_type, color = Source)) +
      geom_point()+
      geom_text_repel() +
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
    
    cor.test(df4$gene_num, df4$IC50, method =  "spearman")
    
    df5 <- df4 %>% filter(Cell_type != "HeLa")
    cor.test(df5$gene_num, df5$IC50, method =  "spearman")
    
    p5 <- df_binded_RNA %>% 
      mutate(Source = ifelse(Cell_type %in% c("PANC-1", "PC-3", "HCT 116", "MG-63", "A2780", "HEK293-T-REx"), "TPM", "RPKM"))  %>%
      ggplot(aes(x = log(RPKM_or_TPM, base = 10), fill = Source)) +
      geom_histogram(bins = 100) +
      theme(plot.title = element_text(face="bold",hjust = 0.5), 
            #legend.position = "none",
            #legend.position = c(0.9, 0.4),
            legend.background = element_rect(color = "black", fill = alpha("white", alpha = 1)),
            panel.grid.major = element_line(colour="gray"),
            panel.grid.minor = element_line(colour="gray", size = 1),
            panel.background = element_blank(), 
            axis.line = element_line(colour="black"),
            axis.text=element_text(size=10,face="bold", color = "black"),
            axis.text.x =element_text(size=10,face="bold", color = "black"),
            axis.text.y =element_text(size=10,face="bold", color = "black"),
            axis.title=element_text(size=10,face="bold", color = "black"),
            aspect.ratio = 0.5
      )+
      facet_wrap(~Cell_type) 
      
  }
}