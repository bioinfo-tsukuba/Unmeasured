pub_chip_plot_CISBP <- function(){
  library(tidyverse)
  TF_information  <- read_tsv("/Users/saeko/MOCCS_paper_public/annotation/CIS-BP/Homo_sapiens_2021_04_21_1_47_am/TF_Information.txt")
  ChIP_annotation <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/experimentList_20230104.tab", col_names = F)
  ChIP_annotation2 <- ChIP_annotation %>% select(X1, X2, X3, X4, X5, X6) %>% filter(X2 == "hg38" & X3 == "TFs and others")
  colnames(ChIP_annotation2) <- c("SRX", "Genome", "Antigen_class", "Antigen", "Cell_type_class", "Cell_type")
  
  # plotし直す
  tgt_DBD_TFs <- TF_information %>% filter(DBDs != "UNKNOWN") %>% .$TF_Name %>% unique()
  
  df3 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/MeSH/Entrez_PMID_AH97924_symbol_human.tsv")
  #Chip_df <- readRDS("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/TF_CTC_sample_num_new.rds")
  Chip_df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/SRX_date_20231004.tsv")
  TF_list <- intersect(unique(Chip_df$Antigen), tgt_DBD_TFs)
  
  tmp  <- df3 %>% filter(Gene_symbol %in% TF_list) %>%
    drop_na(Gene_symbol) %>% select(Gene_symbol, PMID) %>% distinct() %>%
    group_by(Gene_symbol)  %>% summarise(num = n()) 
  
  library(ggrepel)
  tmp2 <- tmp %>% arrange(desc(num))
  top_tf <- tmp2[1:10,]$Gene_symbol %>% as.character()
  tmp3 <- tmp2 %>% mutate(gene_label = ifelse(Gene_symbol %in% top_tf, Gene_symbol, ""))
  
  annotation <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/experimentList.tab", col_names = F)
  annotation2 <- annotation %>% select(X1, X2, X3, X4, X5, X6) %>% distinct()
  colnames(annotation2) <- c("SRX", "genome", "Antigen_class", "Antigen", "Cell_type_class", "Cell_type")
  annotation3 <- annotation2 %>% filter(Antigen %in% tgt_DBD_TFs & genome == "hg38") %>%distinct()
  dim(annotation3)
  
  df_chip <- annotation3 %>% group_by(Antigen) %>% summarise(n_sample = n())
  
  colnames(tmp3) <- c("TF", "pub_num", "gene_label")
  colnames(df_chip) <- c("TF", "chip_num")
  df_join <- df_chip %>% left_join(tmp3, by = "TF") %>% drop_na(pub_num)
  
  cor_test <- cor.test(log(df_join$pub_num, base = 10), log(df_join$chip_num, base = 10), method = "pearson")
  cor <- cor_test$estimate
  p <- df_join %>% ggplot(aes(x = log(pub_num, base = 10), y = log(chip_num,base = 10), label = TF)) +
    geom_point() +
    geom_label_repel()+
    annotate("text", x=1.5, y=2.5, label=cor)+
    xlab("log10(Number of publication)")+
    ylab("log10(Total number of SRX (TF ChIP-seq))")+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_line(colour="gray"),
          panel.grid.minor = element_line(colour="gray", size = 1),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=15,face="bold", color = "black"),
          axis.text.x =element_text(size=15,face="bold", color = "black"),
          axis.text.y =element_text(size=15,face="bold", color = "black"),
          axis.title=element_text(size=15,face="bold", color = "black"),
          aspect.ratio = 1
    )
  
  
  return(p)
  
}