# data check
rm(list=ls())
df_all <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TFmarker/All_TFmarkers.txt")
df_all %>% ggplot(aes(x = 1, fill = `Gene Type`)) +
  geom_bar(stat = "count",
           position = position_stack(), color = "black") +
  geom_text(stat = "count", aes(label = ..count..),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+
  ggtitle("Ratio of TF marker") +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_blank(), 
        axis.text.y =element_blank(), 
        axis.title=element_text(size=12,face="bold"),
        aspect.ratio = 1
  )
View(df_all)

# ChIP-Atlas data
Chip_df1 <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/chip-atlas/experimentList.tab", col_names = F)
Chip_df2 <- Chip_df1 %>%  dplyr::select(X1, X2, X3, X4, X5, X6)
colnames(Chip_df2) <- c("ID", "Genome", "Antigen_class", "Antigen", "Cell_type_class", "Cell_type")
anno_df <- Chip_df2 %>% filter(Genome == "hg38" & Antigen_class == "TFs and others")

TFmarker_ct <- df_all$`Cell Name` %>% unique()
Chip_ct <- anno_df$Cell_type %>% unique()
share_ct <- intersect(TFmarker_ct, Chip_ct)
length(TFmarker_ct)
length(Chip_ct)
length(share_ct)


# Cell Name
df_all$`Cell Name` %>% unique() %>% length()

# Tissue
df_all$`Tissue Type` %>% unique() %>% length()
Tissue <- df_all$`Tissue Type` %>% unique() 
write(Tissue, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TFmarker/All_TFmarkers_tissue.txt")

taiou <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TFmarker/TFmarker_tissue_ChIP-Atlas_Ctc_taiou.tsv")
colnames(taiou) <- c("Tissue_Type", "ChIPAtlas_ctc", "description")
taiou <- taiou %>% select(-description) %>% drop_na(ChIPAtlas_ctc)


# join by Tissue with ChIP-Atlas
colnames(df_all) <- c("PMID", "Gene_Name","Gene_Type","Cell_Name","Cell_Type","Tissue_Type","Experiment_Type","Experimental_Method","Title","Description_of_Gene","Interacting_Gene","CellOntologyID")
df_all2 <- df_all %>% left_join(taiou, by = "Tissue_Type")
df_all2 %>% ggplot(aes(x = 1, fill = ChIPAtlas_ctc)) +
  geom_bar(stat = "count",
           position = position_stack(), color = "black") +
  geom_text(stat = "count", aes(label = ..count..),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+
  ggtitle("Ratio of Cell type class") +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_blank(), 
        axis.text.y =element_blank(), 
        axis.title=element_text(size=12,face="bold"),
        aspect.ratio = 1
  )


# Unmeasured/Measured 
Chip_df <- readRDS("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/TF_CTC_sample_num_new.rds")
TFmarker_tf <- df_all2$Gene_Name %>% unique()
Chip_tf <- anno_df$Antigen %>% unique()
share_tf <- intersect(TFmarker_tf, Chip_tf)
length(TFmarker_tf)
length(Chip_tf)
length(share_tf)
share_ctc <- df_all2 %>% dplyr::filter(ChIPAtlas_ctc != "NA") %>% .$ChIPAtlas_ctc %>% unique()
  
rowname <-rownames(Chip_df)
Chip_df2 <- Chip_df %>% as_tibble() %>% mutate(tf = rowname)
Chip_df_selected <- Chip_df2 %>% select(all_of(share_ctc), tf)
dim(Chip_df_selected)
rowname_selected <- Chip_df_selected$tf
Chip_df_selected2 <- Chip_df_selected %>% select(-tf) %>% as.matrix()
rownames(Chip_df_selected2) <- rowname_selected
dim(Chip_df_selected2)
View(Chip_df_selected2[1:100,])

unmeasured_df <- tibble()
for(j in 1:nrow(Chip_df_selected2)){
  print(paste(j, "/", nrow(Chip_df_selected2)))
  for (i in 1:ncol(Chip_df_selected2)) {
    tgt_tf <- rownames(Chip_df_selected2)[j]
    tgt_ctc <- colnames(Chip_df_selected2)[i]
    tgt_chip_num <- Chip_df_selected2[j, i]
    tgt_row <- df_all2 %>% filter(Gene_Name == tgt_tf & ChIPAtlas_ctc == tgt_ctc) 
    
    if(tgt_chip_num > 1){ #measured
      if(nrow(tgt_row)>0){
        #measured_tf <- c( measured_tf, tgt_row$TF.x)
        #measured_ctc <- c( measured_ctc, tgt_row$Cell_type_class)
      }
    }else{
      if(nrow(tgt_row)>0){  #unmeasured
        print("Unmeasured and TF marker!!")
        if(nrow(unmeasured_df) == 0){
          unmeasured_df  <- tgt_row[,1:13]
        }else{
          unmeasured_df <- unmeasured_df %>% add_row(tgt_row[,1:13])
        }
      }
    }
    rm(tgt_row)
  }
}
write_tsv(unmeasured_df, "/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TFmarker/unmeasured_tf_tissue_TFmarker_selected.tsv")
View(unmeasured_df)
#dim(unmeasured_df)


# TF-markerにはあるのにUnmeasuredだったものの内訳 --------
unmeasured_df <- read_tsv("/Users/saeko/Documents/MOCCS/important_chipseq_prediction/data/TFmarker/unmeasured_tf_tissue_TFmarker_selected.tsv")
unmeasured_df <- unmeasured_df %>% drop_na(ChIPAtlas_ctc)
dim(unmeasured_df)

# 集計用
Chip_tf
tmp <- unmeasured_df %>% filter(Gene_Name %in% Chip_tf) %>% select(Gene_Name, ChIPAtlas_ctc) %>% distinct()
dim(tmp)

head(Chip_df_selected)
tmp2 <- Chip_df_selected %>% pivot_longer(cols = -tf, names_to = "Cell_type_class", values_to = "count") %>% filter(Cell_type_class %in% share_ctc)
tmp3 <- tmp2 %>% filter(count == 0)
dim(tmp3)

# markerごと
## 円グラフ
unmeasured_df %>% ggplot(aes(x = 1, fill = Gene_Type)) +
  geom_bar(stat = "count",
           position = position_stack(), color = "black") +
  geom_text(stat = "count", aes(label = ..count..),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+
  ggtitle("Ratio of TF - Cell type class combinations") +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_blank(), 
        axis.text.y =element_blank(), 
        axis.title=element_text(size=12,face="bold"),
        aspect.ratio = 1
  )

color_list <- rep(c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG")), 30)
unmeasured_df %>% ggplot(aes(x = 1, fill = ChIPAtlas_ctc)) +
  geom_bar(stat = "count",
           position = position_stack(), color = "black") +
  geom_text(stat = "count", aes(label = ..count..),
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values =color_list)+
  coord_polar(theta = "y")+
  ggtitle("Ratio of Cell type") +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_blank(), 
        axis.text.y =element_blank(), 
        axis.title=element_text(size=12,face="bold"),
        aspect.ratio = 1
  )

marker_list <- unique(unmeasured_df$Gene_Type)
p_patch <- c()
for (tgt_marker in marker_list) {
  p <- unmeasured_df %>% filter(Gene_Type == tgt_marker) %>%
    ggplot(aes(x = 1, fill = ChIPAtlas_ctc)) +
    geom_bar(stat = "count",
             position = position_stack(), color = "black") +
    geom_text(stat = "count", aes(label = ..count..),
              position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values =color_list)+
    coord_polar(theta = "y")+
    ggtitle(tgt_marker) +
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_blank(), 
          axis.text.y =element_blank(), 
          axis.title=element_text(size=12,face="bold"),
          aspect.ratio = 1
    )
  if(tgt_marker == marker_list[1]){
    p_patch <- p
  }else{
    p_patch <- p_patch + p
  }
}
plot(p_patch + plot_layout(nrow = 2, ncol = 3), height = 50)


# 棒グラフ
color_list <- rep(c(brewer.pal(10,"Spectral"),brewer.pal(10,"BrBG")), 30)
unmeasured_df %>% 
  ggplot(aes(x = Gene_Type, fill = ChIPAtlas_ctc)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values =color_list)+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=15,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 1
  )

unmeasured_df %>% 
  ggplot(aes(x = Gene_Type, fill = Gene_Name)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values =color_list)+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=15,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 1
  )

unmeasured_df %>% group_by(Gene_Name) %>% summarise(count = n()) %>%
  ggplot(aes(x = reorder(Gene_Name, -count), y = count)) +
  geom_point(size = 1) +
  scale_fill_manual(values =color_list)+
  xlab("Gene name") +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=3,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        aspect.ratio = 0.7
  )

p_patch2 <- c()
for (tgt_marker in marker_list) {
  p <- unmeasured_df %>% group_by(Gene_Name) %>% summarise(count = n()) %>%
    ggplot(aes(x = reorder(Gene_Name, -count), y = count)) +
    geom_point(size = 0.5) +
    ggtitle(tgt_marker) +
    xlab("Gene name") +
    scale_fill_manual(values =color_list)+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=3,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=15,face="bold"),
          axis.title=element_text(size=15,face="bold"),
          aspect.ratio = 0.7
    )
  if(tgt_marker == marker_list[1]){
    p_patch2 <- p
  }else{
    p_patch2 <- p_patch2 + p
  }
}
plot(p_patch2 + plot_layout(nrow = 2, ncol = 3), height = 50)


