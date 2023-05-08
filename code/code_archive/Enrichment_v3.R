library(tidyverse)

# Enrichment解析の結果をプロット ----
SRX_data3 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/ncbi_meta/SRX_date_fix.tsv")
SRX_data3$month <- as.numeric(SRX_data3$month)
SRX_data3$day <- as.numeric(SRX_data3$day)
SRX_data4 <- SRX_data3 %>% arrange(year, month, day, time) %>% mutate(order = 1:nrow(SRX_data3))


# annotation 追加 -----
annotation <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/experimentList.tab", col_names = F)
annotation2 <- annotation %>% select(X1, X2, X3, X4, X5, X6) %>% distinct()
colnames(annotation2) <- c("SRX", "genome", "Antigen_class", "Antigen", "Cell_type_class", "Cell_type")
target_ID_list <- SRX_data3$SRX %>% unique()
length(target_ID_list)
annotation3 <- annotation2 %>% filter(SRX %in% target_ID_list & genome == "hg38") %>%distinct()
dim(annotation3)

SRX_data6 <- SRX_data4 %>% left_join(annotation3, by = "SRX") %>% drop_na(Antigen) %>% 
  arrange(year, month, day, time) 

num_list <- c()
label_list <- c()
for (i in 1:nrow(SRX_data6)) {
  print(paste0(i, "/", nrow(SRX_data6)))
  target_df <- SRX_data6[1:i,]
  tmp_num <- length(unique(target_df$Antigen))
  
  if(i == 1){
    num_list <- tmp_num
  }else{
    num_list <- c(num_list, tmp_num)
  }
  
  # 新規TFかどうかの判定
  tmp <- intersect(unique(target_df[1:(i-1),]$Antigen), target_df[i,]$Antigen)
  if(length(tmp) == 0){
    if(i == 1){
      label_list <- target_df[i,]$Antigen
    }else{
      label_list <- c(label_list, target_df[i,]$Antigen)
    }
  }else{
    if(i == 1){
      label_list <- target_df[i,]$Antigen
    }else{
      label_list <- c(label_list, "")
    }
  }
}
SRX_data7 <- SRX_data6 %>% mutate(TF_uniq_num = num_list, label = label_list)

# Add publication annotation
pub1 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/MeSH/Entrez_PMID_AH97924_symbol_human.tsv")
TF_list <- SRX_data7$Antigen %>% unique()
pub2 <- pub1 %>% filter(Gene_symbol %in% TF_list) %>% select(-descriptorid, -descriptorterm) %>%
  distinct() %>% group_by(Gene_symbol) %>% summarise(pub_count = n())
colnames(pub2) <- c("Antigen", "pub_count")
SRX_data8 <- SRX_data7 %>% left_join(pub2, by = "Antigen")


library(ggside)
SRX_data8 %>% ggplot(aes(x = reorder(Antigen, pub_count), y = factor(year))) +
  geom_point() +
  geom_xsideviolin(aes(y = log(pub_count)), orientation = "y") + 
  geom_xsideboxplot(aes(y = log(pub_count)), orientation = "y") + 
  scale_xsidey_discrete(guide = guide_axis(angle = 0))+
   #geom_ysideviolin(aes(x = pub_count), orientation = "x") +
   #geom_ysideboxplot(aes(x = pub_count), orientation = "x") +
   #scale_ysidex_discrete(guide = guide_axis(angle = 45))+
  theme(plot.title = element_text(face="bold",hjust = 0.5),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=3,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=3,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 1
  )


tgt_TF_list <- SRX_data7 %>% select(SRX, Antigen) %>% 
  filter(Antigen != "Unclassified") %>% filter(Antigen != "GFP") %>% filter(Antigen != "Epitope tags")%>% 
  group_by(Antigen) %>% summarise(n = n()) %>% arrange(desc(n)) %>% .$Antigen
tgt_TF_list <- tgt_TF_list[1:150]

label <- SRX_data7$label %>% unique() 
tmp <- tibble(Antigen = label, order2 = 1:length(label)) %>% filter(label %in% tgt_TF_list) %>% arrange(order2)
tgt_TF_list2 <- tmp$Antigen[1:150]

SRX_data9 <- SRX_data8 %>% filter(Antigen %in% tgt_TF_list2) %>% left_join(tmp, by = "Antigen")

SRX_data9 %>% 
  mutate(order = 1:nrow(SRX_data9)) %>%
  #ggplot(aes(y = reorder(Antigen, pub_count), x = order)) +
  ggplot(aes(y = reorder(Antigen, -order2), x = factor(year), color = pub_count, size = pub_count)) +
  geom_point() +
  scale_size_continuous(range = c(1, 4)) + # サイズの範囲を調整
  xlab("Year") +
  ylab("TF")+
  theme(plot.title = element_text(face="bold",hjust = 0.5),
        #legend.position = "none",
        panel.grid.major = element_line(colour="gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=5,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=5,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 2
  )#+
  #geom_ysideviolin(aes(x = pub_count), orientation = "x", size = 2) +
  #geom_ysideboxplot(aes(x = pub_count), orientation = "x",size = 2) +
  #scale_ysidex_discrete(guide = guide_axis(angle = 45))


label <- SRX_data7$label %>% unique() 
tmp <- tibble(Antigen = label, order2 = 1:length(label)) %>% filter(label %in% tgt_TF_list) %>% arrange(order2)

SRX_data10 <- SRX_data8 %>% left_join(tmp, by = "Antigen")

p <- SRX_data10 %>% 
  mutate(order = 1:nrow(SRX_data10)) %>%
  ggplot(aes(y = reorder(Antigen, pub_count), x = factor(year), color = pub_count, size = log(pub_count))) +
  geom_point() +
  scale_size_continuous(range = c(0, 1)) + # サイズの範囲を調整
  xlab("year") +
  ylab("TF") +
  theme(plot.title = element_text(face="bold",hjust = 0.5),
        #legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=0,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 2
  )+
  #geom_xsidedensity(aes(y = stat(density)), position = "stack") +
  #scale_xsidey_discrete(guide = guide_axis(angle = 0))+
  geom_ysideboxplot(aes(x = pub_count)) +
  scale_ysidex_discrete(guide = guide_axis(angle = 45))

ggsave("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/plot/Enrichment_year_TF.pdf", p)
