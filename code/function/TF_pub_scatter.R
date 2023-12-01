TF_pub_scatter <- function(){
  library(tidyverse)
  df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/MeSH/Entrez_PMID_AH97924_symbol_human.tsv")
  
  ## ChIP-Atlas TF list
  ChIP_annotation <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/experimentList_20230104.tab", col_names = F)
  ChIP_annotation2 <- ChIP_annotation %>% select(X1, X2, X3, X4, X5, X6) %>% filter(X2 == "hg38" & X3 == "TFs and others")
  colnames(ChIP_annotation2) <- c("SRX", "Genome", "Antigen_class", "Antigen", "Cell_type_class", "Cell_type")
  TF_list_ChIP_Atlas <- ChIP_annotation2$Antigen %>% unique()
  df2 <- df %>% filter(Gene_symbol %in% TF_list_ChIP_Atlas) %>%
    drop_na(Gene_symbol) %>% 
    select(Gene_symbol, PMID) %>% 
    distinct() %>%
    group_by(Gene_symbol)  %>% 
    summarise(num = n()) 
  
  library(ggrepel)
  df3 <- df2 %>% arrange(desc(num))
  top_tf <- df3[1:10,]$Gene_symbol %>% as.character()
  df4 <- df3 %>% mutate(gene_label = ifelse(Gene_symbol %in% top_tf, Gene_symbol, ""))
  p1 <- df4 %>%
    ggplot(aes(x = reorder(Gene_symbol, -num), y = num, label = gene_label)) +
    geom_point() +
    geom_text_repel(size = 5, force = 10) +
    xlab("TF") +
    ylab("Number of publications") +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=15,face="bold", colour="black"),
          axis.text.x =element_text(size=0,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=15,face="bold", colour="black"),
          axis.title=element_text(size=15,face="bold", colour="black"),
          aspect.ratio = 0.7
    )
  
  # culumtive plot ----
  count <- 0
  count_list <- c()
  ratio_list <- c()
  for (i in 1:nrow(df3)) {
    tgt_tf <- df3[i,1]
    count <- count + df3[i,2]
    if(i == 1){
      count_list <- count
      ratio_list <- count / sum(df3$num)
    }else{
      count_list <- c(count_list, count) 
      ratio_list <- c(ratio_list, count/sum(df3$num)) 
    }
    
  }
  
  df5 <- df3 %>% mutate(count = as.numeric(count_list), ratio = as.numeric(ratio_list)) %>% mutate(xlabel = 1:nrow(df3))
  top_tf <- df3[1:10,]$Gene_symbol %>% as.character()
  p2 <- df5 %>% mutate(gene_label = ifelse(Gene_symbol %in% top_tf, Gene_symbol, "")) %>%
    #ggplot(aes(reorder(Gene_symbol, -num), y = ratio)) +
    ggplot(aes(x = xlabel, y = ratio, label = gene_label)) +
    geom_point() +
    geom_text_repel(size = 5, force = 7) +
    geom_hline(yintercept=0.25,linetype=2,alpha=0.7,size=1,color='blue') +
    geom_hline(yintercept=0.5,linetype=2,alpha=0.7,size=1,color='red') +
    geom_hline(yintercept=0.75,linetype=2,alpha=0.7,size=1,color='limegreen') +
    geom_vline(xintercept=0.02*nrow(df5),linetype=2,alpha=0.7,size=1,color='blue')+
    geom_vline(xintercept=0.1*nrow(df5),linetype=2,alpha=0.7,size=1,color='red')+
    geom_vline(xintercept=0.25*nrow(df5),linetype=2,alpha=0.7,size=1,color='limegreen')+
    #xlim(c(1, nrow(tmp4)))+
    xlab("TF") +
    ylab("Ratio") +
    annotate("text", x=0.02*nrow(df5),   y= 1, label="2%", size = 8, color = "blue") +
    annotate("text", x=0.1*nrow(df5),   y= 1, label="10%", size = 8, color = "red") +
    annotate("text", x=0.25*nrow(df5),   y= 1, label="25%", size = 8, color = "limegreen") +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=15,face="bold", colour="black"),
          axis.text.x =element_text(size=0,face="bold", colour="black", angle = 45, hjust = 1),
          axis.text.y =element_text(size=15,face="bold", colour="black"),
          axis.title=element_text(size=15,face="bold", colour="black"),
          aspect.ratio = 1
    )
  
  
  
  return(list(p1, p2))
}