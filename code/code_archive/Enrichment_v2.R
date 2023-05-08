# chip-atlas enrichment result
result_enrich <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/all_human_gene_enruchment_result_TFs.tsv", col_names = F)
colnames(result_enrich) <- c("ID", "Antigen class", "TF", "Cell type class", "Cell type", "Number of peaks", "Overlaps/Mydata", "Overlaps/Control",
                             "log_pval", "log_qval", "Fold_enrichment", "FE>1?")
View(result_enrich[1:10,])

res_df2 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/res2.tsv")

# metadata to combine year info
SRA_meta <- read_csv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/ncbi_meta/SraRunInfo_ChIP_human.csv")
head(SRA_meta)

length(unique(SRA_meta$Experiment))
length(unique(result_enrich$ID))

# Release dataを得る ---
SRX_list <- result_enrich$ID %>% unique()
length(SRX_list)
write(SRX_list[1:3], "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/ncbi_meta/test_SRX.txt")
write(SRX_list, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/ncbi_meta/SRX_enrichment_hu38.txt")

# シェルで下記を実行
# esearchとefeatchが使えるようにする：$ conda install -c bioconda entrez-direct
# $cat SRX_enrichment_hu38.txt | xargs -I{} sh -c 'esearch -db sra -query "{}" | efetch -format runinfo' | cut -d',' -f1,2,11 > SRX_date.csv
SRX_date1 <- read_csv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/ncbi_meta/SRX_date.csv", col_names = F)
View(SRX_date1)

# 偶数行のみ抽出
nrow <- nrow(SRX_date1)
SRX_date2 <- SRX_date1 %>% filter(!str_detect(X2, "ReleaseDate")) %>% select(X2, X3) %>% distinct()
colnames(SRX_date2) <- c("Release_date", "SRX")

SRX_date3 <- SRX_date2 %>% separate(Release_date, into = c("year_date", "time"), sep = " ") %>% 
  separate(year_date, into = c("year", "month", "day"), sep = "-") %>%
  select(SRX, year, month, day, time) %>% distinct()
write_tsv(SRX_date3, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/ncbi_meta/SRX_date_fix.tsv")


# Enrichment解析の結果をプロット ----
SRX_data3 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/ncbi_meta/SRX_date_fix.tsv")
SRX_data3$month <- as.numeric(SRX_data3$month)
SRX_data3$day <- as.numeric(SRX_data3$day)

SRX_data4 <- SRX_data3 %>% arrange(year, month, day, time) %>% mutate(order = 1:nrow(SRX_data3))
View(SRX_data4)

# yearごとのSRX数
SRX_data4 %>% group_by(year) %>% summarise(SRX_num =n()) %>% 
  ggplot(aes(x = factor(year), y = SRX_num)) +
  geom_point(size = 5) +
  xlab("Year")+
  ylab("Number of SRX (TF ChIP-seq)")+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_line(colour="gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.7
  )

SRX_data4 %>% group_by(year, month) %>% summarise(SRX_num =n()) %>% 
  unite("year_month", c(year, month), sep = "_") %>%
  ggplot(aes(x = factor(year_month), y = SRX_num)) +
  geom_point(size = 2) +
  xlab("Year_month")+
  ylab("Number of SRX (TF ChIP-seq)")+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_line(colour="gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=6,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.5
  )

# yearごとの累積SRX数
SRX_data5 <- SRX_data4 %>% group_by(year) %>% summarise(SRX_num =n()) 
cum_sum <- cumsum(SRX_data5$SRX_num) #累積和
SRX_data5 <- SRX_data5 %>% mutate(cum = cum_sum)

SRX_data5 %>%
  ggplot(aes(x = factor(year), y = cum)) +
  geom_point(size = 5) +
  xlab("Year")+
  ylab("Total number of SRX (TF ChIP-seq)")+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_line(colour="gray"),
        panel.grid.minor = element_line(colour="gray", size = 1),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.7
  )


  
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

SRX_data7 %>% group_by(year) %>% summarise(TF_uniq_num_max = max(TF_uniq_num)) %>%
  ggplot(aes(x = factor(year), y = TF_uniq_num_max)) +
  geom_point(size = 5) +
  geom_hline(yintercept=1639, linetype="dashed",size=1) +
  xlab("Year")+
  ylab("Total number of unique TFs")+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_line(colour="gray"),
        panel.grid.minor = element_line(colour="gray", size = 1),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.7
  )

library(ggrepel)
SRX_data7 %>% 
  ggplot(aes(x = order, y = TF_uniq_num, label = label)) +
  geom_point(size = 0.5) +
  geom_label_repel(max.overlaps = 20) +
  geom_hline(yintercept=1639, linetype="dashed",size=1) +
  xlab("Time")+
  ylab("Total number of unique TFs")+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_line(colour="gray"),
        panel.grid.minor = element_line(colour="gray", size = 1),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=0,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.7
  )

# ChIP-seqサンプルの多いTF+alphaのみラベル
TF_list_order <- annotation3 %>% filter(Antigen_class == "TFs and others")%>% group_by(Antigen) %>% summarise(n = n()) %>%
  arrange(desc(n)) %>% .$Antigen
top_TF_list <- TF_list_order[1:30]

SRX_data7 %>% group_by(year, month) %>% summarise()
  ggplot(aes(x = order, y = TF_uniq_num, label = label2)) +
  geom_point(size = 0.1) +
  #geom_text() +
  geom_text_repel(max.overlaps = 100) +
  xlab("Time")+
  ylab("Total number of unique TFs")+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_line(colour="gray"),
        panel.grid.minor = element_line(colour="gray", size = 1),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=0,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.7
  )


# TFごとにプロット ----
SRX_data7_2 <- SRX_data7 %>% filter(Antigen %in% top_TF_list) 

SRX_data8 <- c()
for(tgt_tf in top_TF_list){
  print(tgt_tf)
  target_df <- SRX_data7_2 %>% filter(Antigen == tgt_tf) %>% group_by(year) %>% summarise(SRX_num =n()) 
  cum_sum <- cumsum(target_df$SRX_num) #累積和
  target_df2 <- target_df %>% mutate(cum_per_TF = cum_sum, Antigen = tgt_tf)
  
  if(tgt_tf == top_TF_list[1]){
    SRX_data8 <- target_df2
  }else{
    SRX_data8 <- rbind(SRX_data8, target_df2)
  }
}

library(RColorBrewer)
color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG"))
SRX_data8 %>% filter(Antigen != "Epitope tags" & Antigen %in% top_TF_list[1:21]) %>%
  ggplot(aes(x = year, y = cum_per_TF, color = Antigen)) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(values = color_list) +
  xlab("Year")+
  ylab("Total number of SRX (TF ChIP-seq)")+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_line(colour="gray"),
        panel.grid.minor = element_line(colour="gray", size = 1),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.7
  )

# Cell_type_class、SRXとunique TFの累積をプロットする -----
# SRXの累積
SRX_data9 <- c()
for(tgt_ctc in unique(SRX_data7$Cell_type_class)){
  print(tgt_ctc)
  target_df <- SRX_data7 %>% filter(Cell_type_class == tgt_ctc) %>% group_by(year) %>% summarise(SRX_num =n()) 
  cum_sum <- cumsum(target_df$SRX_num) #累積和
  target_df2 <- target_df %>% mutate(cum_per_ctc = cum_sum, Cell_type_class = tgt_ctc)
  
  if(tgt_tf == unique(SRX_data7$Cell_type_class)[1]){
    SRX_data9 <- target_df2
  }else{
    SRX_data9 <- rbind(SRX_data9, target_df2)
  }
}


color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG"))
SRX_data9  %>%
  ggplot(aes(x = year, y = cum_per_ctc, color = Cell_type_class)) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(values = color_list) +
  xlab("Year")+
  ylab("Total number of SRX (TF ChIP-seq)")+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_line(colour="gray"),
        panel.grid.minor = element_line(colour="gray", size = 1),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.7
  )

# TFの累積
SRX_data10 <- c()
for(tgt_ctc in unique(SRX_data7$Cell_type_class)){
  print(tgt_ctc)
  target_df <- SRX_data7 %>% filter(Cell_type_class == tgt_ctc) 
  
  num_list <- c()
  label_list <- c()
  for (i in 1:nrow(target_df)) {
    sub_target_df <- target_df[1:i,]
    tmp_num <- length(unique(sub_target_df$Antigen))
    
    if(i == 1){
      num_list <- tmp_num
    }else{
      num_list <- c(num_list, tmp_num)
    }
    
    # 新規TFかどうかの判定
    tmp <- intersect(unique(sub_target_df[1:(i-1),]$Antigen), sub_target_df[i,]$Antigen)
    if(length(tmp) == 0){
      if(i == 1){
        label_list <- sub_target_df[i,]$Antigen
      }else{
        label_list <- c(label_list, sub_target_df[i,]$Antigen)
      }
    }else{
      if(i == 1){
        label_list <- sub_target_df[i,]$Antigen
      }else{
        label_list <- c(label_list, "")
      }
    }
  }
  target_df2 <- target_df %>% mutate(TF_uniq_num_per_ctc = num_list, label_list_per_ctc = label_list)
  if(tgt_ctc == unique(SRX_data7$Cell_type_class)[1]){
    SRX_data10 <- target_df2
  }else{
    SRX_data10 <- rbind(SRX_data10, target_df2)
  }
}

color_list <- c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG"))
SRX_data10 %>% group_by(year, Cell_type_class) %>% summarise(TF_uniq_num_max = max(TF_uniq_num_per_ctc)) %>%
  ggplot(aes(x = year, y = TF_uniq_num_max, color = Cell_type_class)) +
  geom_point(size = 3) +
  geom_line() +
  xlab("Year")+
  ylab("Total number of unique TFs")+
  scale_color_manual(values = color_list) +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        #legend.position = "none",
        panel.grid.major = element_line(colour="gray"),
        panel.grid.minor = element_line(colour="gray", size = 1),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.7
  )
